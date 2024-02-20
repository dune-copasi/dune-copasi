
#include <dune/copasi/common/tiff_grayscale.hh>

#include <dune/common/exceptions.hh>

#include <gtest/gtest.h>

#include <cassert>
#include <filesystem>

const std::filesystem::path filename_04("data/tiff/flower-minisblack-04.tif");
const std::filesystem::path filename_08("data/tiff/flower-minisblack-08.tif");
const std::filesystem::path filename_16("data/tiff/flower-minisblack-16.tif");

TEST(TestTIFFGrayscale, UnsupportedEncoding)
{
  EXPECT_THROW({ Dune::Copasi::TIFFGrayscale tiff_04(filename_04); }, Dune::NotImplemented);
}

TEST(TestTIFFGrayscale, Compare16vs8Bits)
{
  Dune::Copasi::TIFFGrayscale tiff_08(filename_08);
  Dune::Copasi::TIFFGrayscale tiff_16(filename_16);

  EXPECT_EQ(tiff_08.cols(), tiff_16.cols());
  EXPECT_EQ(tiff_08.rows(), tiff_16.rows());

  short res_unit = 72;
  for (size_t i = 0; i < tiff_08.rows(); i++) {
    for (size_t j = 0; j < tiff_08.cols(); j++) {
      double threshold = 2. / (std::numeric_limits<unsigned char>::max());
      EXPECT_NEAR(tiff_08[i][j] - tiff_16[i][j], 0, threshold);
      double x = (double)j / res_unit;
      double y = (double)i / res_unit;
      EXPECT_NEAR(tiff_08(x, y) - tiff_16(x, y), 0, threshold);
    }
  }
}