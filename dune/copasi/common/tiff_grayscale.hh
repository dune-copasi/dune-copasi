#ifndef DUNE_COPASI_TIFF_GRAYSCALE_HH
#define DUNE_COPASI_TIFF_GRAYSCALE_HH

#include <dune/copasi/common/tiff_file.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <algorithm>
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
    TIFFGrayscaleRow(const TIFFFile& tiff,
                     const T& row)
      : _tiff(tiff)
      , _row(row)
    {
      T* raw_buffer = (T*)_tiff.malloc_scanline();
      auto deleter = [&](auto& ptr) { _tiff.free(ptr); };
      _tiff_buffer = std::shared_ptr<T>(raw_buffer, deleter);
      _tiff.read_scanline(_tiff_buffer.get(), _row);
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
      assert((short)col < _tiff._col_size);
      const T max = std::numeric_limits<T>::max();
      T val = *(_tiff_buffer.get() + col);
      val = _tiff._zero ? val : max - val;
      return (double)val / max;
    }

    /**
     * @brief      The size of this row
     *
     * @return     Number of columns in this row
     */
    std::size_t size() const { return static_cast<std::size_t>(_tiff._col_size); }

    /**
     * @brief      The current row
     *
     * @return     Current row
     */
    std::size_t row() const { return static_cast<std::size_t>(_row); }

  private:
    const TIFFFile& _tiff;
    std::shared_ptr<T> _tiff_buffer;
    const T _row;
  };

public:
  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  filename   The filename
   * @param[in]  max_cache  The maximum row cache.
   */
  TIFFGrayscale(const std::string& filename, std::size_t max_cache = 8)
    : _tiff{filename}
    , _max_cache(max_cache)
  {
    if (_tiff._bits_per_sample != 8 * sizeof(T)) {
      _tiff.close();
      DUNE_THROW(IOError,
                 "TIFF file '" << filename
                               << "' contains a non-readable grayscale field.");
    }
  }

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
    assert((short)row < _tiff._row_size);
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
    int i = static_cast<int>(_tiff._x_res * (x - _tiff._x_off));
    int j = _tiff._row_size - static_cast<int>(_tiff._y_res * (y - _tiff._y_off)) - 1;
    // clamp invalid pixel indices to nearest valid pixel
    i = std::clamp(i, 0, _tiff._col_size - 1);
    j = std::clamp(j, 0, _tiff._row_size - 1);
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
  std::size_t rows() const { return static_cast<std::size_t>(_tiff._row_size); }

  /**
   * @brief      The size of cols for this file
   *
   * @return     Number of cols in this file
   */
  std::size_t cols() const { return static_cast<std::size_t>(_tiff._col_size); }

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
      _row_cache.emplace_back(_tiff, row);

    if (_row_cache.size() >= 8)
      _row_cache.pop_front();

    return *(_row_cache.rbegin());
  }

private:
  TIFFFile _tiff;
  mutable std::deque<TIFFGrayscaleRow> _row_cache;
  const std::size_t _max_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_TIFF_GRAYSCALE_HH
