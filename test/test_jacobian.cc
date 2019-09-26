#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/copasi/gmsh_reader.hh>
#include <dune/copasi/model_diffusion_reaction.hh>
#include <dune/copasi/model_diffusion_reaction.cc>
#include <dune/copasi/model_multidomain_diffusion_reaction.hh>
#include <dune/copasi/enum.hh>

#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

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

    // create a grid
    constexpr int dim = 2;
    using HostGrid = Dune::UGGrid<dim>;
    using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim,1>;
    using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid,MDGTraits>;

    auto& grid_config = config.sub("grid");
    auto level = grid_config.get<int>("initial_level",0);

    log.info("Creating a rectangular grid in {}D"_fmt, dim);

    auto grid_file = grid_config.get<std::string>("file");

    auto [grid_ptr,host_grid_ptr] = Dune::Copasi::GmshReader<Grid>::read(grid_file,config);

    log.debug("Applying global refinement of level: {}"_fmt, level);
    grid_ptr->globalRefine(level);

    auto& model_config = config.sub("model");
    int order = model_config.get<int>("order");
    std::string jacobian = model_config.get<std::string>("jacobian_type");
    
    if (order == 1)
    {
      using Ordering = Dune::PDELab::EntityBlockedOrderingTag;
      constexpr int Order = 1;
      if (jacobian == "analytical")
      {
        constexpr Dune::Copasi::JacobianMethod Jac = Dune::Copasi::JacobianMethod::Analytical;
        using ModelTraits = Dune::Copasi::ModelMultiDomainDiffusionReactionTraits<Grid,Order,Ordering,Jac>;
        Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(grid_ptr,model_config);
        model.run();
      }
      else if (jacobian == "numerical")
      {
        constexpr Dune::Copasi::JacobianMethod Jac = Dune::Copasi::JacobianMethod::Numerical;
        using ModelTraits = Dune::Copasi::ModelMultiDomainDiffusionReactionTraits<Grid,Order,Ordering,Jac>;
        Dune::Copasi::ModelMultiDomainDiffusionReaction<ModelTraits> model(grid_ptr,model_config);
        model.run();
      }
      else 
      {
        DUNE_THROW(Dune::IOError,"Jacobian type " << jacobian << " is not supported by dune-copasi");        
      }
    } 
    else
    {
      DUNE_THROW(Dune::IOError,"Finite element order " << order << " is not supported by dune-copasi");
    }

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
