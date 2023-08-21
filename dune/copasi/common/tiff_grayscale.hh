#ifndef DUNE_COPASI_TIFF_GRAYSCALE_HH
#define DUNE_COPASI_TIFF_GRAYSCALE_HH

#include <dune/copasi/common/tiff_file.hh>

#include <function2/function2.hpp>

#include <functional>
#include <queue>

namespace Dune::Copasi {

/**
 * @brief This class describes a tiff grayscale.
 * @details This class opens, closes, and query information from a single tiff
 * image. Data can be queried using matrix like indices with the bracket
 * operator or using space coordinates using the parethesis operator.
 */
class TIFFGrayscale
{
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
    TIFFGrayscaleRow(const TIFFFile& tiff, const std::size_t& row);
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
    [[nodiscard]] double operator[](const std::size_t& col) const;

    /**
     * @brief      The size of this row
     *
     * @return     Number of columns in this row
     */
    [[nodiscard]] std::size_t size() const;

    /**
     * @brief      The current row
     *
     * @return     Current row
     */
    [[nodiscard]] std::size_t row() const;

  private:
    TIFFFile const* _tiff_ptr;
    std::size_t _row;
    std::size_t _max;
    fu2::unique_function<std::size_t(std::size_t) const noexcept>
      _read_col{};
  };

public:
  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  filename   The filename
   * @param[in]  max_cache  The maximum row cache.
   */
  explicit TIFFGrayscale(const fs::path& filename, std::size_t max_cache = 8);

  /**
   * @brief      Array indexer row operator.
   * @warning    Row object is only a proxy object and ... TODO
   *
   * @param[in]  row   The row index
   *
   * @return     The row
   */
  const TIFFGrayscaleRow& operator[](std::size_t row) const;

  /**
   * @brief      TIFF value by position call operator.
   * @details    Here we assume that TIFFTAG_RESOLUTIONUNIT has same units as
   *             the grid. Since we never check grid units, we also do not check
   *             tiff units
   *
   * @param[in]  pos_x     The x position in same units as the tiff file
   * @param[in]  pos_y     The y position in same units as the tiff file
   *
   *
   * @return     Scaled TIFF value if x and y are in the TIFF domain, 0
   *             otherwise
   */
  [[nodiscard]] double operator()(double pos_x, double pos_y) noexcept;

  /**
   * @brief      The size of rows for this file
   *
   * @return     Number of rows in this file
   */
  [[nodiscard]] std::size_t size() const;

  /**
   * @brief      The size of rows for this file
   *
   * @return     Number of rows in this file
   */
  [[nodiscard]] std::size_t rows() const;

  /**
   * @brief      The size of cols for this file
   *
   * @return     Number of cols in this file
   */
  [[nodiscard]] std::size_t cols() const;

private:
  /**
   * @brief      Cache rows
   *
   * @param[in]  row   The row to be cached
   *
   * @return     The cached row
   */
  const TIFFGrayscaleRow& cache(std::size_t row) const;

  TIFFFile _tiff;
  mutable std::deque<TIFFGrayscaleRow> _row_cache;
  std::size_t _max_cache;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_TIFF_GRAYSCALE_HH
