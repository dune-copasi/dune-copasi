#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/copasi/tiff_grayscale.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <iostream>
#include <set>
#include <cassert>

template<class TGF, class Domain, class Range>
bool test_tiff_grid_function(const TGF& tiff_grid_function, std::map<Domain,Range> points)
{
  return true;
}

int main(int argc, char** argv)
{
  bool failed = false;

  try{
    // initialize mpi
    auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    auto comm = mpi_helper.getCollectiveCommunication();

    // initialize loggers
    // Dune::Logging::Logging::init(comm);

    // using Grid = Dune::YaspGrid<3>;
    // using GridView = typename Grid::LeafGridView;

    // using DF = typename Grid::ctype;
    // using RF = double;
    // int constexpr dimRange = 1;
    // using Range = Dune::FieldVector<RF,dimRange>;

    // Dune::FieldVector<DF,3> L3(1.0); L3[0] = 4;
    // Grid grid(L3,{{4,4}});
    // grid.globalRefine(4);
    // GridView grid_view = grid.leafGridView();
  
    // using TIFFGF = Dune::Copasi::TIFFGridFunction<GridView,RF,dimRange,Range>;


    std::string filename_1("data/tiff/flower-minisblack-04.tif");
    Dune::Copasi::TIFFGrayscale  tiff_1(filename_1); 

    std::string filename_2("data/tiff/flower-minisblack-08.tif");
    Dune::Copasi::TIFFGrayscale  tiff_2(filename_2); 

    for (size_t i = 0; i < 20; i++)
      for (size_t j = 0; j < 20; j++)
        std::cout << (int) tiff_1[j][i] << " " << (int) tiff_2[j][i]/16 << std::endl;
    
    return failed;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
