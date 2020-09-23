#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "grid_function_compare.hh"

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/logging/logging.hh>

#include <ctime>

int
main(int argc, char** argv)
{
  auto stime_c = std::chrono::system_clock::now();
  int end_code = 0;

  // initialize mpi
  auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);

  // Read and parse ini file
  if (argc != 2)
    DUNE_THROW(Dune::IOError, "Wrong number of arguments");
  const std::string config_filename = argv[1];

  Dune::ParameterTree config;
  Dune::ParameterTreeParser::readINITree(config_filename, config);

  // initialize loggers
  auto comm = mpi_helper.getCollectiveCommunication();
  Dune::Logging::Logging::init(comm, config.sub("logging"));
  auto log = Dune::Logging::Logging::logger(config);
  using namespace Dune::Literals;

  try {
    std::time_t stime_t = std::chrono::system_clock::to_time_t(stime_c);
    std::tm stime_tm;
    std::memcpy(&stime_tm, std::localtime(&stime_t), sizeof(std::tm));
    auto stime_s = fmt::format("{:%a %F %T %Z}", stime_tm);
    log.notice("Starting dune-copasi(md) at {}"_fmt, stime_s);
    log.info("Reading configuration file: '{}'"_fmt, config_filename);

    // detailed report of the input paramter tree
    std::stringstream ss;
    config.report(ss);
    log.detail(2, "----"_fmt);
    for (std::string line; std::getline(ss, line);)
      log.detail(2, "{}"_fmt, line);
    log.detail(2, "----"_fmt);

    // create a grid
    constexpr int dim = 2;
    using HostGrid = Dune::UGGrid<dim>;
    using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
    using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

    auto& grid_config = config.sub("grid", true);
    auto level = grid_config.get<int>("initial_level", 0);

    auto grid_file = grid_config.get<std::string>("file");

    auto [md_grid_ptr, host_grid_ptr] =
      Dune::Copasi::MultiDomainGmshReader<Grid>::read(grid_file);

    auto grid_log = Dune::Logging::Logging::componentLogger({}, "grid");
    grid_log.detail("Applying refinement of level: {}"_fmt, level);

    for (int i = 1; i <= level; i++) {
      Dune::Copasi::mark_stripes(*host_grid_ptr);
      md_grid_ptr->preAdapt();
      md_grid_ptr->adapt();
      md_grid_ptr->postAdapt();
    }

    auto& model_config = config.sub("model", true);
    int order = model_config.get<int>("order");

    // create time stepper
    auto timestep_config = model_config.sub("time_stepping", true);
    auto end_time = timestep_config.template get<double>("end");
    auto initial_step = timestep_config.template get<double>("initial_step");
    auto stepper = Dune::Copasi::make_default_stepper(timestep_config);

    auto file = model_config.get("writer.file_path", "");

    // 1. create expression grid functions
    using SubDomainGridView = typename Grid::SubDomainGrid::LeafGridView;
    using GF =
      Dune::Copasi::ExpressionToGridFunctionAdapter<SubDomainGridView, double>;
    auto& compartments_config = model_config.sub("compartments", true);
    const auto& compartments = compartments_config.getValueKeys();
    const std::size_t domains = compartments.size();

    std::vector<std::vector<std::shared_ptr<GF>>> gf_expressions;

    for (std::size_t domain_i = 0; domain_i < domains; ++domain_i) {
      const std::string compartment = compartments[domain_i];
      const std::size_t domain =
        compartments_config.template get<std::size_t>(compartment);
      auto sub_grid_view = md_grid_ptr->subDomain(domain).leafGridView();

      gf_expressions.emplace_back(Dune::Copasi::get_muparser_expressions(
        model_config.sub(compartment + ".compare.expression"), sub_grid_view));
    }

    // get compare paramter tree for a specific variable
    auto setup_param = [&](auto domain, auto var) {
      auto compare_config = model_config.sub(domain + ".compare", true);
      Dune::ParameterTree param;
      // if key is not available for a var, use a fallback on the parent section
      if (compare_config.hasKey("l1_error." + var))
        param["l1_error"] = compare_config["l1_error." + var];
      if (compare_config.hasKey("l2_error." + var))
        param["l2_error"] = compare_config["l2_error." + var];
      if (compare_config.hasKey("linf_error." + var))
        param["linf_error"] = compare_config["linf_error." + var];
      return param;
    };

    auto log_writer = Dune::Logging::Logging::componentLogger({}, "writer");
    std::vector<std::shared_ptr<Dune::VTKSequenceWriter<SubDomainGridView>>>
      writer(domains);
    for (size_t i = 0; i < domains; i++) {
      auto gv = md_grid_ptr->subDomain(i).leafGridView();

      // create directory if necessary
      auto path_entry = fs::directory_entry{ file };
      if (not path_entry.exists()) {
        log_writer.info("Creating output directory '{}'"_fmt,
                        path_entry.path().string());
        std::error_code ec{ 0, std::generic_category() };
        fs::create_directories(path_entry.path(), ec);
        if (ec)
          DUNE_THROW(Dune::IOError,
                    "\n Category: " << ec.category().name() << '\n'
                                    << "Value: " << ec.value() << '\n'
                                    << "Message: " << ec.message() << '\n');
      }

      std::string name = fmt::format(
        "{}-{}-expression", file, compartments_config.getValueKeys()[i]);
      auto base_writer = std::make_shared<Dune::VTKWriter<SubDomainGridView>>(
        gv, Dune::VTK::conforming);
      writer[i] = std::make_shared<Dune::VTKSequenceWriter<SubDomainGridView>>(
        base_writer, name, file, file);
    }

    auto compare_m = [=](const auto& model, const auto& state) {
      // 2. get resulting grid functions
      auto gf_results = model.get_grid_functions(state);
      assert(gf_results.size() == gf_expressions.size());

      // 3. set expression grid functions to current time
      for (auto&& domain_gf_expressions : gf_expressions)
        for (auto&& gf_expression : domain_gf_expressions)
          gf_expression->set_time(state.time);

      // 4. compare expression vs resulting gf
      for (std::size_t domain_i = 0; domain_i < domains; domain_i++) {
        const std::string compartment = compartments[domain_i];
        const std::size_t domain =
          compartments_config.template get<std::size_t>(compartment);
        const auto& diff_vars =
          model_config.sub(compartment + ".diffusion", true).getValueKeys();
        for (std::size_t var_i = 0; var_i < diff_vars.size(); var_i++) {
          auto param_compare = setup_param(compartment, diff_vars[var_i]);
          grid_function_compare(param_compare,
                                *gf_expressions[domain][var_i],
                                *gf_results[domain][var_i]);
          auto vtk_function = Dune::PDELab::makeVTKGridFunctionAdapter(
            gf_expressions[domain][var_i], diff_vars[var_i]);
          writer[domain]->addVertexData(vtk_function);
        }
        if (not file.empty())
          writer[domain]->write(state.time, Dune::VTK::base64);
        writer[domain]->vtkWriter()->clear();
      }

      // write state if requested
      if (not file.empty())
        state.write(file, true);
    };

    if (order == 0) {
      constexpr int Order = 0;
      using ModelTraits =
        Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, Order>;
      Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(
        md_grid_ptr, model_config);
      auto compare = [&](const auto& state) { compare_m(model, state); };
      compare(model.state()); // compare initial condition
      stepper.evolve(model, initial_step, end_time, compare);
    } else if (order == 1) {
      constexpr int Order = 1;
      using ModelTraits =
        Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, Order>;
      Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(
        md_grid_ptr, model_config);
      auto compare = [&](const auto& state) { compare_m(model, state); };
      compare(model.state()); // compare initial condition
      stepper.evolve(model, initial_step, end_time, compare);
    } else if (order == 2) {
      constexpr int Order = 2;
      using ModelTraits =
        Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, Order>;
      Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(
        md_grid_ptr, model_config);
      auto compare = [&](const auto& state) { compare_m(model, state); };
      compare(model.state()); // compare initial condition
      stepper.evolve(model, initial_step, end_time, compare);
    } else {
      DUNE_THROW(Dune::IOError,
                 "Finite element order " << order
                                         << " is not supported by dune-copasi");
    }
  } catch (Dune::Exception& e) {
    log.error("Dune reported error:"_fmt);
    log.error(2, "{}"_fmt, e.what());
    end_code = 1;
  } catch (std::exception& e) {
    log.error("C++ reported error:"_fmt);
    log.error(2, "{}"_fmt, e.what());
    end_code = 1;
  } catch (...) {
    log.error("Unknown exception thrown!"_fmt);
    end_code = 1;
  }

  if (end_code)
    log.notice("dune-copasi(md) finished with some errors :("_fmt);
  else
    log.notice("dune-copasi(md) successfully finished :)"_fmt);
  return end_code;
}
