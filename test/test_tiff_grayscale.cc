#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/common/tiff_grayscale.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <cassert>

int
main(int argc, char** argv)
{
  bool failed = false;

  try {

    std::string filename_04("data/tiff/flower-minisblack-04.tif");
    std::string filename_08("data/tiff/flower-minisblack-08.tif");
    std::string filename_16("data/tiff/flower-minisblack-16.tif");

    static_assert(sizeof(unsigned char) == 1);
    static_assert(sizeof(unsigned short) == 2);
    Dune::Copasi::TIFFGrayscale<unsigned char> tiff_08(filename_08);
    Dune::Copasi::TIFFGrayscale<unsigned short> tiff_16(filename_16);

    try {
      Dune::Copasi::TIFFGrayscale<unsigned char> tiff_08(filename_04);
      failed |= true;
    } catch (...) {
    }

    try {
      Dune::Copasi::TIFFGrayscale<unsigned char> tiff_08(filename_16);
      failed |= true;
    } catch (...) {
    }

    try {
      Dune::Copasi::TIFFGrayscale<unsigned short> tiff_16(filename_04);
      failed |= true;
    } catch (...) {
    }

    try {
      Dune::Copasi::TIFFGrayscale<unsigned short> tiff_16(filename_08);
      failed |= true;
    } catch (...) {
    }

    assert(tiff_08.cols() == tiff_16.cols());
    assert(tiff_08.rows() == tiff_16.rows());

    short res_unit = 72;

    for (size_t i = 0; i < tiff_08.rows(); i++) {
      for (size_t j = 0; j < tiff_08.cols(); j++) {
        double threshold = 1. / (std::numeric_limits<unsigned char>::max());
        failed |= abs(tiff_08[i][j] - tiff_16[i][j]) > threshold;
        double x = (double)j / res_unit;
        double y = (double)i / res_unit;
        failed |= abs(tiff_08(x, y) - tiff_16(x, y)) > threshold;
      }
    }

    return failed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
