#ifdef HAVE_DUNE_COPASI_CONFIG_H
#include <dune/copasi/config.h>
#endif

#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <ctime>
#include <vector>

static std::string string_help =
  "Usage: dune-copasi-md CONFIG_FILE\n"
  "\n"
  "Execute numerical simulation of reaction-diffusion systems on multiple\n"
  "compartments. The CONFIG_FILE is a DUNE INI file with the\n"
  "parameters to perform the simulation.\n";

int
main(int argc, char** argv)
{
  std::vector<std::string> cmd_line_args(argv, argv+argc);

  if (cmd_line_args.size() != 2) {
    std::cerr << string_help;
    return 1;
  }

  if (cmd_line_args[1] == "--help" or cmd_line_args[1] == "-h") {
    std::cout << string_help;
    return 0;
  }

  auto stime_c = std::chrono::system_clock::now();
  int end_code = 0;

  // initialize mpi
  auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);

  // Read and parse ini file
  const std::string config_filename = cmd_line_args[1];
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
    auto evolve_model = [](const auto& config, auto dim, auto order) {
      using namespace Dune::Copasi;
      using Dune::TypeTree::treePath;

      using HostGrid = Dune::UGGrid<dim>;
      using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
      using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

      auto& grid_config = config.sub("grid", true);
      auto& model_config = config.sub("model", true);
      auto& timestep_config = model_config.sub("time_stepping", true);

      auto grid_file = grid_config.template get<std::string>("file");
      auto level = grid_config.template get<int>("initial_level");

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

      auto model = DiffusionReactionModel{model_config};
      auto state = model.make_multidomain_multicomponent_state(md_grid_ptr, order);
      state.time = timestep_config.template get<double>("begin");

      for (std::size_t domain = 0; domain < state.space().degree(); ++domain){
        // for each sub-domain interpret initial value expressions into a grid function...
        auto& domain_space = state.space().child(domain);
        auto domain_name = domain_space.name();
        auto& initial = model_config.sub(domain_name).sub("initial", true);
        auto mc_gf = make_multicomponent_grid_function(initial, domain_space.child(0).gridView(), false);
        add_tiff_to_grid_function(mc_gf, model_config.sub("data"));
        mc_gf.setTime(state.time);
        // ...and interpolate it to the current state of the sub-domain
        model.interpolate(state, mc_gf, treePath(domain));
      }

      {
        std::vector writers{state.writer};
        writers.resize(state.space().degree());
        for (std::size_t domain = 0; domain < state.space().degree(); ++domain)
          writers[domain] = model.make_vtk_writer(treePath(domain));

        state.writer = [=](const auto& state, const fs::path& path, bool append) mutable {
          for (std::size_t domain = 0; domain < state.space().degree(); ++domain) {
            auto domain_name = state.space().child(domain).name();
            auto domain_path = path.string() + "-" + domain_name;
            writers[domain](state, domain_path, append);
          }
        };
      }

      // setup writer
      std::string file = model_config.get("writer.file_path", "");
      auto at_end_of_step = [=](const auto& state) {
        if (not file.empty())
          state.write(file, true);
      };
      at_end_of_step(state); // write current time

      // create time stepper
      auto end_time = timestep_config.template get<double>("end");
      auto initial_step = timestep_config.template get<double>("initial_step");
      auto stepper = make_default_stepper(timestep_config);

      // auto state_out = state;
      stepper.evolve(model, state, state, initial_step, end_time, at_end_of_step);
    };

    // maximum polynomial order instantiated
    constexpr auto max_order = 2;
    if (order == 0 or order > max_order)
      DUNE_THROW(Dune::IOError,
                 "Finite element order " << order
                                         << " is not supported by dune-copasi");

    // static loop that instantiates 2D and 3D polynomial orders but runs the
    // dynamic ones
    auto order_range = Dune::range(Dune::Indices::_1, Dune::index_constant<max_order + 1>{});
    Dune::Hybrid::forEach(order_range, [&](auto static_order) {
      // instantiate models with 2D grids
      if (dim == 2 and order == static_order)
        evolve_model(config, Dune::Indices::_2, static_order);

      // instantiate models with 3D grids
      if (dim == 3 and order == static_order)
#ifdef DUNE_COPASI_COMPILE_3D
        evolve_model(config, Dune::Indices::_3, static_order);
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
