#ifndef DUNE_COPASI_TIFF_FILE_HH
#define DUNE_COPASI_TIFF_FILE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <tiffio.h>

#include <memory>
#include <string>

namespace Dune::Copasi {

/**
 * @brief Simple tiff file interface
 */
struct TIFFFile
{
  //! opens a tiff file from the file system
  TIFFFile(const std::string& filename);

  //! Destructor
  ~TIFFFile();

  //! Closes the file
  void close();

  //! Allocates buffer to scan lines of the file
  void* malloc_scanline() const;

  //! Deallocates tiff type pointers
  void free(void* ptr) const;

  //! Reads one line on the buffer for a given row
  void read_scanline(void* ptr, std::size_t row) const;

  void* _tiff_file;
  bool _zero;
  short _bits_per_sample;
  short _row_size;
  short _col_size;
  float _x_res, _x_off;
  float _y_res, _y_off;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_TIFF_FILE_HH
