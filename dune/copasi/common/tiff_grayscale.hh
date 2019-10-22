#ifndef DUNE_COPASI_TIFF_GRAYSCALE_HH
#define DUNE_COPASI_TIFF_GRAYSCALE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <tiffio.h>

#include <memory>
#include <queue>
#include <string>

namespace Dune::Copasi {

/**
 * @brief      This class describes a tiff grayscale.
 * @details    This class opens, closes, and quere information from a single
 *             tiff image. Data can be querid using matrix like indices with the
 *             bracket operator or using space coordinates using the parethesis
 *             operator.
 *
 * @tparam     T     Unsigned integer type to represent each pixel information.
 */
template<class T>
class TIFFGrayscale
{
  static_assert(std::is_integral_v<T>, "T must be an integral type");
  static_assert(std::is_unsigned_v<T>, "T must be an unsigned type");

  /**
   * @brief      This class represents a single row of a tiff grayscale
   */
  struct TIFFGrayscaleRow
  {
    /**
     * @brief      Constructs a tiff grayscale row
     *
     * @param      tiff_file  The tiff file
     * @param[in]  row        The row
     * @param[in]  col_size   The column size
     * @param[in]  zero       Zero based grayscale?
     */
    TIFFGrayscaleRow(TIFF* const tiff_file,
                     const T& row,
                     const short& col_size,
                     const bool& zero)
      : _row(row)
      , _col_size(col_size)
      , _zero(zero)
    {
      T* raw_buffer = (T*)_TIFFmalloc(TIFFScanlineSize(tiff_file));
      auto deleter = [](auto& ptr) { _TIFFfree(ptr); };
      _tiff_buffer = std::shared_ptr<T>(raw_buffer, deleter);
      TIFFReadScanline(tiff_file, _tiff_buffer.get(), _row);
    }

    /**
     * @brief      Array indexer column operator.
     * @details    This operator always scales the gayscale value with its
     *             maximum value. This way, the results is independent of the
     *             template argument T being used
     *
     * @param[in]  col   The column for the current row
     *
     * @return     The result of the array indexer scaled between 0 and 1.
     */
    double operator[](const T& col) const
    {
      assert((short)col < _col_size);
      const T max = std::numeric_limits<T>::max();
      T val = *(_tiff_buffer.get() + col);
      val = _zero ? val : max - val;
      return (double)val / max;
    }

    /**
     * @brief      The size of this row
     *
     * @return     Number of columns in this row
     */
    std::size_t size() const { return static_cast<std::size_t>(_col_size); }

    /**
     * @brief      The current row
     *
     * @return     Current row
     */
    std::size_t row() const { return static_cast<std::size_t>(_row); }

  private:
    std::shared_ptr<T> _tiff_buffer;
    const T _row;
    const short _col_size;
    const bool _zero;
  };

public:
  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  filename   The filename
   * @param[in]  max_cache  The maximum row cache.
   */
  TIFFGrayscale(const std::string& filename, std::size_t max_cache = 8)
    : _tiff_file(TIFFOpen(filename.c_str(), "r"))
    , _max_cache(max_cache)
  {
    if (not _tiff_file)
      DUNE_THROW(IOError, "Error opening TIFF file '" << filename << "'.");

    short photometric;
    TIFFGetField(_tiff_file, TIFFTAG_PHOTOMETRIC, &photometric);
    if ((photometric != PHOTOMETRIC_MINISWHITE) and
        (photometric != PHOTOMETRIC_MINISBLACK))
      DUNE_THROW(IOError,
                 "TIFF file '" << filename << "' must be in grayscale.");
    _zero = (bool)photometric;

    short bits_per_sample;
    TIFFGetField(_tiff_file, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);

    if (bits_per_sample != 8 * sizeof(T)) {
      TIFFClose(_tiff_file);
      DUNE_THROW(IOError,
                 "TIFF file '" << filename
                               << "' contains a non-readable grayscale field.");
    }

    TIFFGetField(_tiff_file, TIFFTAG_IMAGELENGTH, &_row_size);
    TIFFGetField(_tiff_file, TIFFTAG_IMAGEWIDTH, &_col_size);
    TIFFGetField(_tiff_file, TIFFTAG_XRESOLUTION, &_x_res);
    TIFFGetField(_tiff_file, TIFFTAG_YRESOLUTION, &_y_res);
    assert(FloatCmp::gt(_x_res, 0.f));
    assert(FloatCmp::gt(_y_res, 0.f));
    _x_off = _y_off = 0.;
    TIFFGetField(_tiff_file, TIFFTAG_XPOSITION, &_x_off);
    TIFFGetField(_tiff_file, TIFFTAG_YPOSITION, &_y_off);
  }

  /**
   * @brief      Destroys the object.
   */
  ~TIFFGrayscale() { TIFFClose(_tiff_file); }

  /**
   * @brief      Array indexer row operator.
   * @warning    If rows are read concurrently, this object has to be copied
   *
   * @param[in]  row   The row index
   *
   * @return     The row
   */
  const TIFFGrayscaleRow& operator[](T row) const
  {
    assert((short)row < _row_size);
    return cache(row);
  }

  /**
   * @brief      TIFF value by position call operator.
   * @details    Here we assume that TIFFTAG_RESOLUTIONUNIT has same units as
   *             the grid. Since we never check grid units, we also do not check
   *             tiff units
   *
   * @param[in]  x     The x position in same units as the tiff file
   * @param[in]  y     The x position in same units as the tiff file
   *
   * @tparam     DF    Domain Field type
   *
   * @return     Scaled TIFF value if x and y are in the TIFF domain, 0
   *             otherwise
   */
  template<class DF>
  double operator()(const DF& x, const DF& y)
  {
    // return 0 if not in the domain
    if (FloatCmp::lt((float)x, _x_off) or FloatCmp::lt((float)y, _y_off))
      return 0;

    const T i = _x_res * (_x_off + x);
    const T j = _row_size - _y_res * (_y_off + y) - 1;

    // return 0 if not in the domain
    if (i >= cols() or j >= rows())
      return 0;
    else
      return (*this)[j][i];
  }

  /**
   * @brief      The size of rows for this file
   *
   * @return     Number of rows in this file
   */
  std::size_t size() const { return rows(); }

  /**
   * @brief      The size of rows for this file
   *
   * @return     Number of rows in this file
   */
  std::size_t rows() const { return static_cast<std::size_t>(_row_size); }

  /**
   * @brief      The size of cols for this file
   *
   * @return     Number of cols in this file
   */
  std::size_t cols() const { return static_cast<std::size_t>(_col_size); }

private:
  /**
   * @brief      Cache rows
   *
   * @param[in]  row   The row to be cached
   *
   * @return     The cached row
   */
  const TIFFGrayscaleRow& cache(T row) const
  {
    auto it = _row_cache.rbegin();
    while ((it != _row_cache.rend()) and (it->row() != row))
      it++;

    if (it != _row_cache.rend())
      return *it;
    else
      _row_cache.emplace_back(_tiff_file, row, _col_size, _zero);

    if (_row_cache.size() >= 8)
      _row_cache.pop_front();

    return *(_row_cache.rbegin());
  }

private:
  TIFF* _tiff_file;
  mutable std::deque<TIFFGrayscaleRow> _row_cache;
  short _row_size;
  short _col_size;
  float _x_res, _x_off;
  float _y_res, _y_off;
  bool _zero;
  const std::size_t _max_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_TIFF_GRAYSCALE_HH