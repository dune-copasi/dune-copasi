#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/filesystem.hh>
#include <dune/copasi/common/tiff_file.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <tiff.h>
#include <tiffio.h>

#include <cstddef>
#include <cstdint>

namespace Dune::Copasi {

TIFFFile::TIFFFile(const fs::path& filename)
  : _tiff_file(static_cast<TIFF*>(TIFFOpen(filename.c_str(), "r")))
{
  TIFF* tiff_file = static_cast<TIFF*>(_tiff_file);
  if (tiff_file == nullptr) {
    throw format_exception(IOError{}, "Error opening TIFF file '{}'", filename.string());
  }

  uint16_t photometric = 0;
  // NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,hicpp-vararg): this is how the library is used
  TIFFGetField(tiff_file, TIFFTAG_PHOTOMETRIC, &photometric);
  if ((photometric != PHOTOMETRIC_MINISWHITE) and (photometric != PHOTOMETRIC_MINISBLACK)) {
    throw format_exception(IOError{}, "TIFF file '{}' must be in grayscale", filename.string());
  }
  _info.zero = static_cast<bool>(photometric);

  TIFFGetField(tiff_file, TIFFTAG_BITSPERSAMPLE, &_info.bits_per_sample);
  TIFFGetField(tiff_file, TIFFTAG_IMAGELENGTH, &_info.row_size);
  TIFFGetField(tiff_file, TIFFTAG_IMAGEWIDTH, &_info.col_size);
  TIFFGetField(tiff_file, TIFFTAG_XRESOLUTION, &_info.x_res);
  TIFFGetField(tiff_file, TIFFTAG_YRESOLUTION, &_info.y_res);
  if (FloatCmp::le(_info.x_res, 0.F) or FloatCmp::le(_info.y_res, 0.F)) {
    throw format_exception(IOError{}, "TIFF file '{}' has negative resolution", filename.string());
  }
  TIFFGetField(tiff_file, TIFFTAG_XPOSITION, &_info.x_off);
  TIFFGetField(tiff_file, TIFFTAG_YPOSITION, &_info.y_off);
  // NOLINTEND(cppcoreguidelines-pro-type-vararg,hicpp-vararg)
}

TIFFFile::~TIFFFile()
{
  if (_tiff_file != nullptr) {
    close();
  }
}

void
TIFFFile::close()
{
  TIFFClose(static_cast<TIFF*>(_tiff_file));
  _tiff_file = nullptr;
}

auto
TIFFFile::malloc_scanline() const -> void*
{
  return _TIFFmalloc(TIFFScanlineSize(static_cast<TIFF*>(_tiff_file)));
}

void
TIFFFile::free(void* ptr)
{
  _TIFFfree(ptr);
}

void
TIFFFile::read_scanline(void* ptr, std::size_t row) const
{
  TIFFReadScanline(static_cast<TIFF*>(_tiff_file), ptr, row);
}

} // namespace Dune::Copasi
