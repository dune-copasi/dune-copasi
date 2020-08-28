#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/grid/mark_stripes.hh>
#include <dune/copasi/grid/multidomain_gmsh_reader.hh>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

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
    Dune::ParameterTreeParser::readINITree(config_filename, config);

    // initialize loggers
    Dune::Logging::Logging::init(comm, config.sub("logging"));

    using namespace Dune::Literals;
    auto log = Dune::Logging::Logging::logger(config);
    log.notice("Starting dune-copasi"_fmt);

    log.debug("Input config file '{}':"_fmt, config_filename);

    Dune::Logging::LoggingStream ls(false, log.indented(2));
    if (log.level() >= Dune::Logging::LogLevel::debug)
      config.report(ls);

    // Grid setup
    constexpr int dim = 2;
    using HostGrid = Dune::UGGrid<dim>;
    using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
    using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

    auto& grid_config = config.sub("grid");
    auto level = grid_config.get<int>("initial_level", 0);

    auto grid_file = grid_config.get<std::string>("file");

    auto [md_grid_ptr, host_grid_ptr] =
      Dune::Copasi::MultiDomainGmshReader<MDGrid>::read(grid_file, config);

    using Grid = typename MDGrid::SubDomainGrid;
    using GridView = typename Grid::Traits::LeafGridView;

    log.debug("Applying refinement of level: {}"_fmt, level);

    for (int i = 1; i <= level; i++) {
      Dune::Copasi::mark_stripes(*host_grid_ptr);

      md_grid_ptr->preAdapt();
      md_grid_ptr->adapt();
      md_grid_ptr->postAdapt();
    }

    auto& model_config = config.sub("model");
    auto& compartments_map = model_config.sub("compartments");
    int order = model_config.get<int>("order");

    // TODO: Use OS for different domains when is ready
    if (compartments_map.getValueKeys().size() != 1)
      DUNE_THROW(
        Dune::NotImplemented,
        "Multiple compartments per model are not allowed in this executable");

    // get subdomain grid as a shared pointer
    auto compartment = compartments_map.getValueKeys().front();
    int domain = compartments_map.template get<int>(compartment);
    std::shared_ptr<Grid> grid_ptr =
      Dune::stackobject_to_shared_ptr(md_grid_ptr->subDomain(domain));

    if (order == 0) {
      constexpr int Order = 0;
      using ModelTraits =
        Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, Order>;
      Dune::Copasi::ModelDiffusionReaction<ModelTraits> model(grid_ptr,
                                                              model_config);
      model.run();
    } else if (order == 1) {
      constexpr int Order = 1;
      using ModelTraits =
        Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, Order>;
      Dune::Copasi::ModelDiffusionReaction<ModelTraits> model(grid_ptr,
                                                              model_config);
      model.run();
    } else if (order == 2) {
      constexpr int Order = 2;
      using ModelTraits =
        Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, Order>;
      Dune::Copasi::ModelDiffusionReaction<ModelTraits> model(grid_ptr,
                                                              model_config);
      model.run();
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
