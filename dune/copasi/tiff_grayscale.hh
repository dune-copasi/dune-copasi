#ifndef DUNE_COPASI_TIFF_GRAYSCALE_HH
#define DUNE_COPASI_TIFF_GRAYSCALE_HH

#include <dune/common/exceptions.hh>

#include <tiffio.h>

#include <string>

namespace Dune::Copasi {

class TIFFGrayscale
{

  struct TIFFGrayscaleRow
  {
    TIFFGrayscaleRow(const unsigned char* const tiff_buffer, const bool& half_word)
      : _tiff_buffer(tiff_buffer)
      , _half_word(half_word)
    {}

    inline unsigned char operator[] (std::size_t col) const
    {
      if (not _half_word)
        return *(_tiff_buffer+col);
      else if (col%2 == 0)
        return (*(_tiff_buffer+col/2) & 0b1111'0000) >> 4;
      else
        return (*(_tiff_buffer+(col/2) )) & 0b0000'1111;
    }

    const unsigned char* const _tiff_buffer;
    const bool& _half_word;
  };

public:

  TIFFGrayscale(const std::string& filename)
    : _tiff_file(TIFFOpen(filename.c_str(), "r"))
  {
    if (not _tiff_file)
      DUNE_THROW(IOError,"Error opening TIFF file '" << filename << "'.");

    short photometric;
    TIFFGetField(_tiff_file,TIFFTAG_PHOTOMETRIC,&photometric);
    if ((photometric != PHOTOMETRIC_MINISBLACK) and (photometric != PHOTOMETRIC_MINISBLACK))
      DUNE_THROW(IOError,"TIFF file '" << filename << "' must be in grayscale.");


    short bits_per_sample;
    TIFFGetField(_tiff_file,TIFFTAG_BITSPERSAMPLE,&bits_per_sample);
    if ((bits_per_sample != 4) and (bits_per_sample != 8))
      DUNE_THROW(IOError,"TIFF file '" << filename << "' contains a non-readable grayscale field.");

    _half_word = (bits_per_sample==4);

    _tiff_buffer = (unsigned char*) _TIFFmalloc(TIFFScanlineSize(_tiff_file));
  }

  inline TIFFGrayscaleRow operator[] (std::size_t row) const
  {
    TIFFReadScanline(_tiff_file, _tiff_buffer, row);
    return TIFFGrayscaleRow(_tiff_buffer,_half_word);
  }

  ~TIFFGrayscale()
  {
    _TIFFfree(_tiff_buffer);
    TIFFClose(_tiff_file);
  }

private:
  TIFF* _tiff_file;
  unsigned char* _tiff_buffer;
  bool _half_word;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_TIFF_GRAYSCALE_HH