#ifndef DUNE_COPASI_GRID_CELL_DATA_HH
#define DUNE_COPASI_GRID_CELL_DATA_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/grid/concepts/gridview.hh>

#include <optional>
#include <string>
#include <vector>

namespace Dune::Copasi {

/**
 * @brief Container for cell data of a grid view
 *
 * @tparam GV   The grid view type
 * @tparam T    The data type to store
 */
template<Dune::Concept::GridView GV, class T>
class CellData
{
  // Index mapper for the grid view entities
  using IndexMapper = MultipleCodimMultipleGeomTypeMapper<GV>;
  using IndexSet = typename GV::IndexSet;

  // Allow makeLeafGridViewCellData to see the contents of the level grid view case 
  friend CellData<typename GV::Grid::LevelGridView, T>;

public:
  using GridView = GV;

  // Entity
  using Entity = typename GridView::template Codim<0>::Entity;

  explicit CellData(GridView grid_view)
    : _index_mapper{ grid_view, mcmgElementLayout() }
  {
    // extract the level of the grid view if any
    if (std::same_as<GridView, typename GridView::Grid::LevelGridView>) {
      auto entity_it = grid_view.template begin<0>();
      if (entity_it != grid_view.template end<0>())
        _grid_view_level.emplace(0);
    }
  }

  //! creates a cell data object with a copy data for the leaf grid view
  /**
   * @brief Create cell data for leaf grid view
   * @details Makes a copy of this CellData into another CellData that operates in the
   * leaf grid view
   * 
   * @return std::unique_ptr<CellData<typename GV::Grid::LeafGridView, T>> Cell data for leaf grid view
   */
  std::unique_ptr<CellData<typename GV::Grid::LeafGridView, T>> makeLeafGridViewCellData() const {
    using LeafGridView = typename GV::Grid::LeafGridView;
    if constexpr (std::same_as<GV, LeafGridView>) {
      return std::make_unique<CellData<LeafGridView, T>>(*this);
    } else {
      if (_grid_view_level and _grid_view_level.value() != 0)
        throw format_exception(
          InvalidStateException{},
          "Leaf grid view cell data can only be made out of other leaf grid view (copy) or a 0-th level grid view (projection)");
      
      // copy values level grid view values into leaf grid view values
      auto leaf_gv = _index_mapper.gridView().grid().leafGridView();
      auto leaf_cd = std::make_unique<CellData<LeafGridView, T>>(leaf_gv);
      leaf_cd->reserve(_keys.size());
      auto lvl_cap = capacity();
      auto leaf_cap = leaf_cd->capacity();
      auto sz = size();
      for (const auto& entity : elements(leaf_gv)) {
        auto lvl_offset = index(entity) * lvl_cap;
        auto leaf_offset = leaf_cd->index(entity) * leaf_cap;
        for (std::size_t id = 0; id != sz; ++id) {
          bool mask = leaf_cd->_cell_mask[leaf_offset + id] = _cell_mask[lvl_offset + id];
          leaf_cd->_cell_data[leaf_offset + id] = mask ? _cell_data[lvl_offset + id] : std::numeric_limits<T>::quiet_NaN();
        }
      }
      leaf_cd->_keys = _keys;
      return leaf_cd;
    }
  }

  //! Number of elements per entity that the container has currently allocated space for
  std::size_t capacity() const noexcept { return _cell_data.size() / _index_mapper.size(); }

  //! Number of elements per entity in the container
  std::size_t size() const noexcept { return _keys.size(); }

  //! Number of elements per entity that the container can ever hold container
  std::size_t max_size() const noexcept { return _cell_data.max_size() / _index_mapper.size(); }

  //! The size of the gridview - This gives the number of elements in the grid view
  std::size_t gridview_size() const noexcept { return _index_mapper.size(); }

  //! Return the keys of the cell data
  const std::vector<std::string>& keys() const noexcept { return _keys; }

  //! Reservers new_cap elements per entity in the container
  // note that this does not honor the standard guarantees meaning that if an exception is thrown,
  // there is no guarantee that the data will be unaffected
  void reserve(std::size_t new_cap)
  {
    if (capacity() < new_cap)
      realloc(new_cap);
  }

  void shrink_to_fit()
  {
    if (auto sz = size(); sz != capacity())
      realloc(sz);
  }

  // contiguous range
  template<std::ranges::input_range R,
           class Proj = std::identity,
           std::indirect_unary_predicate<std::projected<std::ranges::iterator_t<R>, Proj>> Pred>
    requires std::assignable_from<T&, std::ranges::range_value_t<R>>
  void addData(std::string_view key, R&& r, Pred pred, Proj proj = {})
  {
    auto ei = size();
    grow(size() + 1);
    pushKey(key);
    auto cap = capacity();
    auto it = r.begin();
    for (std::size_t i = 0; i != _index_mapper.size(); ++i, ++it) {
      bool mask = _cell_mask[i * cap + ei] = std::invoke(pred, std::invoke(proj, *it));
      _cell_data[i * cap + ei] = mask ? *it : std::numeric_limits<T>::quiet_NaN();
    }
    if (it != r.end())
      throw format_exception(
        InvalidStateException{},
        "Range has a size of '{}', but there are '{}' entities in the container",
        std::ranges::size(r),
        _index_mapper.size());
  }

  template<std::ranges::input_range R>
  requires std::assignable_from<T&, std::ranges::range_value_t<R>>
  void addData(std::string_view key, R&& r)
  {
    addData(
      key, std::forward<R>(r), [](auto&&) { return std::true_type{}; }, std::identity{});
  }

  // map range
  template<std::ranges::input_range R,
           class Proj = std::identity,
           std::indirect_unary_predicate<std::projected<std::ranges::iterator_t<R>, Proj>> Pred>
    requires std::assignable_from<T&, std::tuple_element_t<1, std::ranges::range_value_t<R>>>
  void addData(std::string_view key, R&& r, Pred pred, Proj proj = {})
  {
    auto ei = size();
    grow(size() + 1);
    pushKey(key);
    auto cap = capacity();
    auto imsz = _index_mapper.size();
    for (const auto& [i, val] : r) {
      if(i >= _index_mapper.size())
        throw format_exception(RangeError{}, "Key '{}' assinges a values in an index '{}' bigger than the size of the grid view '{}'", key, i, imsz);
      bool mask = _cell_mask[i * cap + ei] = std::invoke(pred, std::invoke(proj, val));
      _cell_data[i * cap + ei] = mask ? val : std::numeric_limits<T>::quiet_NaN();
    }
  }

  // map range
  template<std::ranges::input_range R>
  requires std::assignable_from<T&, std::tuple_element_t<1, std::ranges::range_value_t<R>>>
  void addData(std::string_view key, R&& r)
  {
    addData(
      key, std::forward<R>(r), [](auto&&) { return std::true_type{}; }, std::identity{});
  }

  // copies data and data mask into buffers
  void getData(const Concept::IndexableEntity<IndexSet> auto& entity,
               std::vector<T>& cell_data,
               std::vector<bool>& cell_mask) const
  {
    auto entity_index = index(entity);
    auto cap = capacity();
    auto sz = size();
    cell_data.resize(sz);
    cell_mask.resize(sz);
    auto offset = entity_index * cap;
    for (std::size_t id = 0; id != sz; ++id) {
      bool mask = cell_mask[id] = _cell_mask[offset + id];
      cell_data[id] = mask ? _cell_data[offset + id] : std::numeric_limits<T>::quiet_NaN();
    }
  }

  // id is the index with respect to keys
  std::optional<T> getData(std::size_t id,
                           const Concept::IndexableEntity<IndexSet> auto& entity) const
  {
    auto entity_index = index(entity);
    auto data_id = entity_index * capacity() + id;
    return _cell_mask[data_id] ? std::make_optional<T>(_cell_data[data_id]) : std::nullopt;
  }

  // slower version
  std::optional<T> getData(std::string_view key,
                           const Concept::IndexableEntity<IndexSet> auto& entity) const
  {
    auto key_it = std::find(_keys.begin(), _keys.end(), key);
    if (key_it == _keys.end())
      throw format_exception(RangeError{}, "Key '{}' does not exist in cell data manager", key);
    auto id = std::distance(_keys.begin(), key_it);
    return getData(id, entity);
  }

private:
  template<Concept::IndexableEntity<IndexSet> E>
  auto index(const E& entity) const
  {
    if (_grid_view_level and (*_grid_view_level != entity.level())) {
      E _entity = entity;
      while (_entity.level() > *_grid_view_level) {
        _entity = _entity.father();
      }
      if (_entity.level() != *_grid_view_level)
        throw format_exception(InvalidStateException{},
                               "Invoking a cell data manager with an entity of lower level than "
                               "the grid view attached to it is an invalid operation!");
      return _index_mapper.index(_entity);
    } else {
      return _index_mapper.index(entity);
    }
  }

  void pushKey(std::string_view key)
  {
    auto key_it = std::find(_keys.begin(), _keys.end(), key);
    if (key_it != _keys.end())
      throw format_exception(
        InvalidStateException{}, "Key '{}' is already present in the cell data", key);
    _keys.push_back(std::string{ key });
  }

  void grow(std::size_t new_size)
  {
    if (new_size < capacity())
      return; // no-op
    if (auto ms = max_size(); new_size > ms / 2)
      reserve(ms); // maximum possible size
    else
      reserve(std::max(new_size, size() * 2)); // textbook growth factor of 2!
  }

  void realloc(std::size_t new_cap)
  {
    auto cap = capacity();
    if (new_cap < size())
      throw std::length_error{ "Requested raw capacity is smaller than current size" };
    if (new_cap > max_size())
      throw std::length_error{ "Requested capacity is bigger than maximum size" };

    std::vector<T> new_cell_data;
    new_cell_data.reserve(new_cap * _index_mapper.size());
    auto data_begin = _cell_data.begin();

    std::vector<bool> new_cell_mask;
    new_cell_mask.reserve(new_cap * _index_mapper.size());
    auto mask_begin = _cell_mask.begin();

    for (std::size_t i = 0; i != _index_mapper.size(); ++i) {
      // move data
      auto data_end = std::next(data_begin, cap);
      new_cell_data.insert(new_cell_data.end(),
                           std::make_move_iterator(data_begin),
                           std::make_move_iterator(data_end));
      new_cell_data.resize(new_cap * (i + 1), std::numeric_limits<T>::quiet_NaN());
      data_begin = data_end;
      // move mask
      auto mask_end = std::next(mask_begin, cap);
      new_cell_mask.insert(new_cell_mask.end(),
                           std::make_move_iterator(mask_begin),
                           std::make_move_iterator(mask_end));
      new_cell_mask.resize(new_cap * (i + 1), false);
      mask_begin = mask_end;
    }
    std::swap(_cell_data, new_cell_data);
    std::swap(_cell_mask, new_cell_mask);
  }

  IndexMapper _index_mapper;

  // the level of the grid view, if none, leaf grid view or empty level grid view
  std::optional<int> _grid_view_level;

  std::vector<std::string> _keys;

  std::vector<T> _cell_data;
  std::vector<bool> _cell_mask;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_CELL_DATA_HH
