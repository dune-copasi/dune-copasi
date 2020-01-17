#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "grid_function_compare.hh"

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>
#include <dune/logging/loggingstream.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

template<class Model>
void
model_compare(const Dune::ParameterTree& param, Model& model)
{
  auto do_step = [&]() {
    return Dune::FloatCmp::lt(model.current_time(), model.end_time());
  };

  // 1. create expression grid functions
  auto grid = model.states().begin()->second.grid;
  using Grid = std::decay_t<decltype(*grid)>;
  static_assert(Dune::Copasi::Concept::isMultiDomainGrid<Grid>());
  using SubDomainGridView = typename Grid::SubDomainGrid::LeafGridView;
  using GF =
    Dune::Copasi::ExpressionToGridFunctionAdapter<SubDomainGridView, double>;
  const auto& compartments = param.sub("compartments").getValueKeys();
  const std::size_t domains = compartments.size();

  std::vector<std::vector<std::shared_ptr<GF>>> gf_expressions;

  for (std::size_t domain_i = 0; domain_i < domains; ++domain_i) {
    const std::string compartement = compartments[domain_i];
    const std::size_t domain =
      param.sub("compartments").template get<std::size_t>(compartement);
    auto sub_grid_view = grid->subDomain(domain).leafGridView();

    gf_expressions.emplace_back(Dune::Copasi::get_muparser_expressions(
      param.sub(compartments[domain] + ".compare.expression"), sub_grid_view));
  }

  auto compare = [&]() {
    auto setup_param = [&](auto domain, auto var) {
      auto compare_config = param.sub(domain + ".compare");
      Dune::ParameterTree param;

      if (compare_config.hasKey("l1_error." + var))
        param["l1_error"] = compare_config["l1_error." + var];
      if (compare_config.hasKey("l2_error." + var))
        param["l2_error"] = compare_config["l2_error." + var];
      if (compare_config.hasKey("linf_error." + var))
        param["linf_error"] = compare_config["linf_error." + var];
      return param;
    };

    // 2. get resulting grid functions
    auto gf_results = model.get_grid_functions();
    assert(gf_results.size() == gf_expressions.size());

    // 3. set expression grid functions to current time
    for (auto&& domain_gf_expressions : gf_expressions)
      for (auto&& gf_expression : domain_gf_expressions)
        gf_expression->set_time(model.current_time());

    // 4. compare expression vs resulting gf
    for (std::size_t domain_i = 0; domain_i < domains; domain_i++) {
      const std::string compartement = compartments[domain_i];
      const std::size_t domain =
        param.sub("compartments").template get<std::size_t>(compartement);
      const auto& diff_vars =
        param.sub(compartement + ".diffusion").getValueKeys();
      for (std::size_t var_i = 0; var_i < diff_vars.size(); var_i++) {
        auto param_compare = setup_param(
          param.sub("compartments").getValueKeys()[domain_i], diff_vars[var_i]);
        grid_function_compare(param_compare,
                              *gf_expressions[domain][var_i],
                              *gf_results[domain][var_i]);
      }
    }
  };

  compare();

  while (do_step()) {
    model.step();

    compare();

    if (model.adaptivity_policy() != Dune::Copasi::AdaptivityPolicy::None)
      if (do_step()) {
        model.mark_grid();
        model.pre_adapt_grid();
        model.adapt_grid();
        model.post_adapt_grid();
      }
  }
}

int
main(int argc, char** argv)
{

  try {
    // initialize mpi
    auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    auto comm = mpi_helper.getCollectiveCommunication();

    // Read and parse ini file
    if (argc != 2)
      DUNE_THROW(Dune::IOError, "Wrong number of arguments");
    const std::string config_filename = argv[1];

    Dune::ParameterTree config;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree(config_filename, config);

    // initialize loggers
    Dune::Logging::Logging::init(comm, config.sub("logging"));

    using namespace Dune::Literals;
    auto log = Dune::Logging::Logging::logger(config);
    log.notice("Starting dune-copasi"_fmt);

    log.debug("Input config file '{}':"_fmt, config_filename);

    Dune::Logging::LoggingStream ls(false, log.indented(2));
    if (log.level() >= Dune::Logging::LogLevel::debug)
      config.report(ls);

    // create a grid
    constexpr int dim = 2;
    using HostGrid = Dune::UGGrid<dim>;
    using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
    using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

    auto& grid_config = config.sub("grid");
    auto level = grid_config.get<int>("initial_level", 0);

    auto grid_file = grid_config.get<std::string>("file");

    auto [md_grid_ptr, host_grid_ptr] =
      Dune::Copasi::MultiDomainGmshReader<Grid>::read(grid_file, config);

    log.debug("Applying refinement of level: {}"_fmt, level);

    for (int i = 1; i <= level; i++) {
      Dune::Copasi::mark_stripes(*host_grid_ptr);

      md_grid_ptr->preAdapt();
      md_grid_ptr->adapt();
      md_grid_ptr->postAdapt();
    }

    auto& model_config = config.sub("model");
    int order = model_config.get<int>("order");

    if (order == 0) {
      constexpr int Order = 0;
      using ModelTraits =
        Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, Order>;
      Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(
        md_grid_ptr, model_config);
      model_compare(model_config, model);
    } else if (order == 1) {
      constexpr int Order = 1;
      using ModelTraits =
        Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, Order>;
      Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(
        md_grid_ptr, model_config);
      model_compare(model_config, model);
    } else if (order == 2) {
      constexpr int Order = 2;
      using ModelTraits =
        Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, Order>;
      Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(
        md_grid_ptr, model_config);
      model_compare(model_config, model);
    } else {
      DUNE_THROW(Dune::IOError,
                 "Finite element order " << order
                                         << " is not supported by dune-copasi");
    }
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
  return 0;
}
