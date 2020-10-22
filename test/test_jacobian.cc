#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.cc>
#include <dune/copasi/model/diffusion_reaction.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.cc>
#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <ctime>
#include <random>

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

    int order = config.template get<int>("model.order");
    int dim = config.template get<int>("grid.dimension");

    if (dim != 2 and dim != 3)
      DUNE_THROW(Dune::IOError, "Only 2D and 3D grids are alloed!");

    // lambda that instantiates and evolves a multidomain model
    auto compare_jacobian = [&log](const auto& config, auto dim, auto order) {
      using namespace Dune::Copasi;

      using HostGrid = Dune::UGGrid<dim>;
      using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
      using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

      auto& grid_config = config.sub("grid", true);
      auto level = grid_config.template get<int>("initial_level");

      auto grid_file = grid_config.template get<std::string>("file");

      auto [md_grid_ptr, host_grid_ptr] =
        MultiDomainGmshReader<Grid>::read(grid_file);

      auto grid_log = Dune::Logging::Logging::componentLogger({}, "grid");
      grid_log.detail("Applying refinement of level: {}"_fmt, level);

      for (int i = 1; i <= level; i++) {
        mark_stripes(*host_grid_ptr);
        md_grid_ptr->preAdapt();
        md_grid_ptr->adapt();
        md_grid_ptr->postAdapt();
      }

      auto& model_config = config.sub("model", true);

      using AnalyticTraits = ModelMultiDomainP0PkDiffusionReactionTraits<
        Grid,
        order,
        Dune::PDELab::EntityBlockedOrderingTag,
        JacobianMethod::Analytical>;
      ModelMultiDomainDiffusionReaction<AnalyticTraits> analytic(md_grid_ptr,
                                                                 model_config);

      using NumericTraits = ModelMultiDomainP0PkDiffusionReactionTraits<
        Grid,
        order,
        Dune::PDELab::EntityBlockedOrderingTag,
        JacobianMethod::Numerical>;
      ModelMultiDomainDiffusionReaction<NumericTraits> numeric(md_grid_ptr,
                                                               model_config);

      auto analytic_grid_op = analytic.get_stationary_grid_operator();
      auto numeric_grid_op = numeric.get_stationary_grid_operator();

      using AnalyticGOP = std::decay_t<decltype(*analytic_grid_op)>;
      using NumericGOP = std::decay_t<decltype(*numeric_grid_op)>;

      using AnalyticMatrix = typename AnalyticGOP::Traits::Jacobian;
      using NumericMatrix = typename NumericGOP::Traits::Jacobian;

      auto analytic_matrix = AnalyticMatrix{*analytic_grid_op};
      auto numeric_matrix = NumericMatrix{*numeric_grid_op};

      auto compare_matrix = [&](const auto& an_x, const auto& nu_x) {
        auto log_matrix = [log](const auto& name, const auto& mat) {
          if (log.level() >= Dune::Logging::LogLevel::trace) {
            std::stringstream ss;
            Dune::printmatrix(ss,mat,name,"");
            for (std::string line; std::getline(ss, line);)
              log.trace(2, "{}"_fmt, line);
          }
        };

        using namespace Dune::PDELab::Backend;
        log.info("Comparing jacobians"_fmt);
        log.info(2, "Vector analytic norm: {}"_fmt,native(an_x).two_norm());
        log.info(2, "Vector numeric norm: {}"_fmt,native(nu_x).two_norm());
        numeric_matrix = 0.;
        analytic_matrix = 0.;

        analytic_grid_op->jacobian(an_x,analytic_matrix);
        numeric_grid_op->jacobian(nu_x,numeric_matrix);

        auto analitic_norm = native(analytic_matrix).frobenius_norm();
        auto numeric_norm = native(numeric_matrix).frobenius_norm();
        log.info(2, "Analytic jacobian norm (an_norm): {}"_fmt,analitic_norm);
        log.info(2, "Numeric jacobian norm (num_norm): {}"_fmt,numeric_norm);

        log_matrix("Analytic Matrix", native(analytic_matrix));
        log_matrix("Numeric Matrix", native(numeric_matrix));

        static_assert(std::is_same_v<AnalyticMatrix, NumericMatrix>);
        native(analytic_matrix) -= native(numeric_matrix);
        log_matrix("Diff Matrix", native(analytic_matrix));
        auto diff_norm = native(analytic_matrix).frobenius_norm();

        // max allowed relative jacobian error
        auto max_error = model_config.get("jacobian_error",1e-7);

        log.info(2, "Difference norm (diff_norm): {}"_fmt,diff_norm);
        log.info(2, "Error (diff_norm/num_norm): {}"_fmt,diff_norm/numeric_norm);
        log.info(2, "Max allowed error: {}"_fmt,max_error);
        if (Dune::FloatCmpOps<double>{}.gt(diff_norm/numeric_norm, max_error))
        {
          log.error("Max allowed error exceded!"_fmt);
          DUNE_THROW(Dune::MathError, "Max allowed error exceded!");
        }
      };

      auto& analytic_x = *analytic.state().coefficients;
      auto& numeric_x = *numeric.state().coefficients;

      // Make the two backends point to the same vector
      static_assert(std::is_same_v<decltype(analytic_x), decltype(numeric_x)>);
      analytic_x.attach(numeric_x.storage());

      analytic_x = 0.;
      compare_matrix(analytic_x,numeric_x);

      analytic_x = 1.;
      compare_matrix(analytic_x,numeric_x);

      {
        std::uniform_real_distribution<double> unif(0.,10.);
        std::default_random_engine re;
        auto rand = [&](){return unif(re);};
        std::generate(analytic_x.begin(), analytic_x.end(), rand);
        compare_matrix(analytic_x,numeric_x);
      }

      {
        std::uniform_real_distribution<double> unif(0.,10000.);
        std::default_random_engine re;
        auto rand = [&](){return unif(re);};
        std::generate(analytic_x.begin(), analytic_x.end(), rand);
        compare_matrix(analytic_x,numeric_x);
      }
    };

    // maximum polynomial order instantiated
    constexpr auto max_order = 2;
    if (order > max_order)
      DUNE_THROW(Dune::IOError,
                 "Finite element order " << order
                                         << " is not supported by dune-copasi");

    // static loop that instantiates 2D and 3D polynomial orders but runs the
    // dynamic ones
    auto order_range = Dune::range(Dune::index_constant<max_order + 1>{});
    Dune::Hybrid::forEach(order_range, [&](auto static_order) {
      // instantiate models with 2D grids
      if (dim == 2 and order == static_order)
        compare_jacobian(config, Dune::Indices::_2, static_order);

      // instantiate models with 3D grids
      if (dim == 3 and order == static_order)
#ifdef DUNE_COPASI_COMPILE_3D
        compare_jacobian(config, Dune::Indices::_3, static_order);
#else
        DUNE_THROW(Dune::IOError,"dune_copasi was not built with 3D support");
#endif
    });

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
