#ifndef DUNE_COPASI_TIFF_FILE_HH
#define DUNE_COPASI_TIFF_FILE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <filesystem>
#include <memory>

namespace Dune::Copasi {

// NOLINTBEGIN(altera-struct-pack-align)

/**
 * @brief Simple tiff file interface
 */
struct TIFFFile
{

  //! Basic information about the tiff file
  struct Info
  {
    uint32_t row_size = 0;
    uint32_t col_size = 0;
    float x_res = 0., x_off = 0.;
    float y_res = 0., y_off = 0.;
    uint16_t bits_per_sample = 0;
    bool zero = false;
  };

  //! opens a tiff file from the file system
  explicit TIFFFile(const std::filesystem::path& filename);

  TIFFFile(const TIFFFile&) = delete;
  TIFFFile(TIFFFile&&) = delete;

  TIFFFile& operator=(const TIFFFile&) = delete;
  TIFFFile& operator=(TIFFFile&&) = delete;

  //! Destructor, closes the file
  ~TIFFFile();

  //! Closes the file
  void close();

  //! Allocates buffer to scan lines of the file
  [[nodiscard]] void* malloc_scanline() const;

  //! Deallocates tiff type pointers
  static void free(void* ptr);

  //! Reads one line on the buffer for a given row
  void read_scanline(void* ptr, std::size_t row) const;

  [[nodiscard]] const Info& info() const { return _info; }

private:
  void* _tiff_file;
  Info _info;
};

// NOLINTEND(altera-struct-pack-align)

} // namespace Dune::Copasi

#endif // DUNE_COPASI_TIFF_FILE_HH
