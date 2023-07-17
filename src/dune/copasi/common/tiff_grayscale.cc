#include <dune/copasi/common/tiff_grayscale.hh>

#include <dune/copasi/common/tiff_file.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <fmt/core.h>

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>

namespace Dune::Copasi {

TIFFGrayscale::TIFFGrayscaleRow::TIFFGrayscaleRow(const TIFFFile& tiff, const std::size_t& row)
  : _tiff_ptr(&tiff)
  , _row(row)
  , _max{ std::size_t{ 1 } << _tiff_ptr->info().bits_per_sample }
{
  auto* raw_buffer = static_cast<std::byte*>(_tiff_ptr->malloc_scanline());
  auto tiff_buffer =
    std::unique_ptr<std::byte, decltype(&TIFFFile::free)>(raw_buffer, &TIFFFile::free);
  _tiff_ptr->read_scanline(tiff_buffer.get(), _row);
  std::size_t bytes = _tiff_ptr->info().bits_per_sample / CHAR_BIT;

  // depending on the encoding, we need to interpret ponters differently
  using Types = std::tuple<uint8_t, uint16_t, uint32_t, uint64_t>;
  Hybrid::forEach(std::make_index_sequence<std::tuple_size_v<Types>>{}, [&](auto type_id) {
    using Type = std::tuple_element_t<type_id, Types>;
    if (bytes == sizeof(Type) and tiff_buffer) {
      _read_col = [bytes,
                   _tiff_buffer = std::move(tiff_buffer)](std::size_t col) noexcept -> std::size_t {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return *reinterpret_cast<Type* const>(_tiff_buffer.get() + (bytes * col));
      };
    }
  });
  if (tiff_buffer) {
    DUNE_THROW(
      NotImplemented,
      fmt::format("Encoding with {} bits not implemented", _tiff_ptr->info().bits_per_sample));
  }
}

auto
TIFFGrayscale::TIFFGrayscaleRow::operator[](const std::size_t& col) const -> double
{
  assert(col < _tiff_ptr->info().col_size);
  std::size_t val = _read_col(col);
  val = _tiff_ptr->info().zero ? val : _max - val;
  return static_cast<double>(val) / static_cast<double>(_max);
}

auto
TIFFGrayscale::TIFFGrayscaleRow::size() const -> std::size_t
{
  return _tiff_ptr->info().col_size;
}

auto
TIFFGrayscale::TIFFGrayscaleRow::row() const -> std::size_t
{
  return _row;
}

TIFFGrayscale::TIFFGrayscale(const fs::path& filename, std::size_t max_cache)
  : _tiff{ filename }
  , _max_cache(max_cache)
{
}

auto
TIFFGrayscale::operator[](std::size_t row) const -> const TIFFGrayscale::TIFFGrayscaleRow&
{
  assert(row < _tiff.info().row_size);
  return cache(row);
}

[[nodiscard]] auto
TIFFGrayscale::operator()(double pos_x, double pos_y) noexcept -> double
{
  auto row = static_cast<uint32_t>(_tiff.info().x_res * (pos_x - _tiff.info().x_off));
  auto col = _tiff.info().row_size -
             static_cast<uint32_t>(_tiff.info().y_res * (pos_y - _tiff.info().y_off)) - 1;
  // clamp invalid pixel indices to nearest valid pixel
  row = std::clamp(row, uint32_t{ 0 }, _tiff.info().col_size - 1);
  col = std::clamp(col, uint32_t{ 0 }, _tiff.info().row_size - 1);
  return (*this)[col][row];
}

[[nodiscard]] auto
TIFFGrayscale::size() const -> std::size_t
{
  return rows();
}

[[nodiscard]] auto
TIFFGrayscale::rows() const -> std::size_t
{
  return _tiff.info().row_size;
}

[[nodiscard]] auto
TIFFGrayscale::cols() const -> std::size_t
{
  return _tiff.info().col_size;
}

auto
TIFFGrayscale::cache(std::size_t row) const -> const TIFFGrayscale::TIFFGrayscaleRow&
{
  auto cahce_it = _row_cache.rbegin();
  while ((cahce_it != _row_cache.rend()) and (cahce_it->row() != row)) {
    cahce_it++;
  }

  if (cahce_it != _row_cache.rend()) {
    return *cahce_it;
  }

  _row_cache.emplace_back(_tiff, row);
  if (_row_cache.size() >= _max_cache) {
    _row_cache.pop_front();
  }

  return *(_row_cache.rbegin());
}

} // namespace Dune::Copasi
