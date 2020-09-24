#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/common/stepper.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

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
    log.notice("Starting dune-copasi(sd) at {}"_fmt, stime_s);
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

      using HostGrid = Dune::UGGrid<dim>;
      using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
      using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

      auto& grid_config = config.sub("grid", true);
      auto level = grid_config.template get<int>("initial_level");
      auto grid_file = grid_config.template get<std::string>("file");

      auto [md_grid_ptr, host_grid_ptr] =
        MultiDomainGmshReader<MDGrid>::read(grid_file);

      auto grid_log = Dune::Logging::Logging::componentLogger({}, "grid");
      grid_log.detail("Applying refinement of level: {}"_fmt, level);

      for (int i = 1; i <= level; i++) {
        mark_stripes(*host_grid_ptr);
        md_grid_ptr->preAdapt();
        md_grid_ptr->adapt();
        md_grid_ptr->postAdapt();
      }

      auto& model_config = config.sub("model", true);
      // TODO: Use OS for different domains when is ready
      auto& compartments_map = model_config.sub("compartments", true);
      if (compartments_map.getValueKeys().size() != 1)
        DUNE_THROW(
          Dune::NotImplemented,
          "Multiple compartments per model are not allowed in this executable");

      // get subdomain grid as a shared pointer
      auto compartment = compartments_map.getValueKeys().front();
      int domain = compartments_map.template get<int>(compartment);
      auto grid_ptr =
        Dune::stackobject_to_shared_ptr(md_grid_ptr->subDomain(domain));

      using Grid = typename MDGrid::SubDomainGrid;
      using GridView = typename Grid::Traits::LeafGridView;
      using Traits = ModelP0PkDiffusionReactionTraits<Grid, GridView, order>;
      ModelDiffusionReaction<Traits> model(grid_ptr, model_config);

      // setup writer
      auto file = model_config.get("writer.file_path", "");
      auto write_output = [=](const auto& state) {
        if (not file.empty())
          state.write(file, true);
      };
      write_output(model.state()); // write initial condition

      // create time stepper
      auto timestep_config = model_config.sub("time_stepping", true);
      auto end_time = timestep_config.template get<double>("end");
      auto initial_step = timestep_config.template get<double>("initial_step");
      auto stepper = Dune::Copasi::make_default_stepper(timestep_config);

      stepper.evolve(model, initial_step, end_time, write_output);
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
        evolve_model(config, Dune::Indices::_2, static_order);

      // instantiate models with 3D grids
      if (dim == 3 and order == static_order)
        evolve_model(config, Dune::Indices::_3, static_order);
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
    log.notice("dune-copasi(sd) finished with some errors :("_fmt);
  else
    log.notice("dune-copasi(sd) successfully finished :)"_fmt);
  return end_code;
}
