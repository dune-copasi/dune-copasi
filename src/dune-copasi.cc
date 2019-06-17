#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

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
    const std::string inifilename = argv[1];

    Dune::ParameterTree inifile;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree(inifilename, inifile);


    // initialize loggers
    Dune::Logging::Logging::init(comm,inifile.sub("logging"));

    using namespace Dune::Literals;
    auto log = Dune::Logging::Logging::logger();
    log.notice("Starting dune-copasi"_fmt);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
