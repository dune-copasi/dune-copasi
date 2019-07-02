#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/copasi/model_diffusion_reaction.hh>
#include <dune/copasi/model_diffusion_reaction.cc>

#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <iostream>

int main(int argc, char** argv)
{

  try{
    // initialize mpi
    auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    auto comm = mpi_helper.getCollectiveCommunication();

    // Read and parse ini file
    if (argc!=2)
      DUNE_THROW(Dune::IOError, "Wrong number of arguments");
    const std::string config_filename = argv[1];

    Dune::ParameterTree config;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree(config_filename, config);

    // initialize loggers
    Dune::Logging::Logging::init(comm,config.sub("logging"));

    using namespace Dune::Literals;
    auto log = Dune::Logging::Logging::logger(config);
    log.notice("Starting dune-copasi"_fmt);


    // test the current code... 

    // create a grid
    constexpr int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    using Domain = Dune::FieldVector<double,2>;

    auto& grid_config = config.sub("grid");
    auto level = grid_config.get<int>("initial_level",0);
    auto upper_right = grid_config.get<Domain>("extensions",{1.,1.});
    auto elements = grid_config.get<std::array<uint, 2>>("cells",{10,10});

    log.info("Creating a rectangular grid in {}D"_fmt, dim);
    log.debug("Grid extensions: {}"_fmt, upper_right);
    log.debug("Grid cells: {}"_fmt, elements);

    Domain origin(.0);

    std::shared_ptr<Grid> grid;
    grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(origin,
                                                             upper_right,
                                                             elements);

    log.debug("Applying global refinement of level: {}"_fmt, level);
    grid->globalRefine(level);

    // instantiate a model
    auto& model_config = config.sub("model");
    Dune::Copasi::ModelDiffusionReaction model(grid,model_config);

    model.run();

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
