#include <dune/copasi/common/tiff_file.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <tiffio.h>

#include <string>

namespace Dune::Copasi {

TIFFFile::TIFFFile(const std::string& filename)
  : _tiff_file((TIFF*)TIFFOpen(filename.c_str(), "r"))
{
  TIFF* tiff_file = (TIFF*)_tiff_file;
  if (not tiff_file)
    DUNE_THROW(IOError, "Error opening TIFF file '" << filename << "'.");

  short photometric;
  TIFFGetField(tiff_file, TIFFTAG_PHOTOMETRIC, &photometric);
  if ((photometric != PHOTOMETRIC_MINISWHITE) and
      (photometric != PHOTOMETRIC_MINISBLACK))
    DUNE_THROW(IOError, "TIFF file '" << filename << "' must be in grayscale.");
  _zero = (bool)photometric;

  TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &_bits_per_sample);
  TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &_row_size);
  TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &_col_size);
  TIFFGetField(tiff_file, TIFFTAG_XRESOLUTION, &_x_res);
  TIFFGetField(tiff_file, TIFFTAG_YRESOLUTION, &_y_res);
  assert(FloatCmp::gt(_x_res, 0.f));
  assert(FloatCmp::gt(_y_res, 0.f));
  _x_off = _y_off = 0.;
  TIFFGetField(tiff_file, TIFFTAG_XPOSITION, &_x_off);
  TIFFGetField(tiff_file, TIFFTAG_YPOSITION, &_y_off);
}

TIFFFile::~TIFFFile()
{
  if (_tiff_file)
    close();
}

void
TIFFFile::close()
{
  TIFFClose((TIFF*)_tiff_file);
  _tiff_file = nullptr;
}

void*
TIFFFile::malloc_scanline() const
{
  return _TIFFmalloc(TIFFScanlineSize((TIFF*)_tiff_file));
}

void
TIFFFile::free(void* ptr) const
{
  _TIFFfree(ptr);
}

void
TIFFFile::read_scanline(void* ptr, std::size_t row) const
{
  TIFFReadScanline((TIFF*)_tiff_file, ptr, row);
}

} // namespace Dune::Copasi
