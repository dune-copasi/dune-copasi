#ifndef DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH
#define DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH

#include <dune-copasi-config.h>

#include <dune/grid/common/capabilities.hh>

#include <dune/geometry/type.hh>

#include <dune/common/exceptions.hh>

namespace Dune::Copasi {

/**
 * @brief      Determines if grid view has a single geometry type.
 * @details    If the grid only knows about one geometry type at compile-time,
 *             this function is constexpr and always true. Notice however,
 *             that the countrary does not hold and the check is done at run-time.
 *
 * @param[in]  grid_view  The grid view
 *
 * @tparam     GridView   The grid view type
 *
 * @return     True if single geometry type, False otherwise.
 */
template<class GridView>
std::enable_if_t<Capabilities::hasSingleGeometryType<typename GridView::Grid>::v,
                 bool> constexpr has_single_geometry_type(const GridView& /*unused*/)
{
  return true;
}

/**
 * @brief      Determines if grid view has a single geometry type.
 * @details    If the grid only knows about one geometry type at compile-time,
 *             this function is constexpr and always true. Notice however,
 *             that the countrary does not hold and the check is done at run-time.
 *
 * @param[in]  grid_view  The grid view
 *
 * @tparam     GridView   The grid view type
 *
 * @return     True if single geometry type, False otherwise.
 */
template<class GridView>
std::enable_if_t<not Capabilities::hasSingleGeometryType<typename GridView::Grid>::v,bool>
has_single_geometry_type(const GridView& grid_view)
{
  if (grid_view.size(0) == 0) {
    return true;
  }
  GeometryType gt = grid_view.template begin<0>()->geometry().type();
  for (auto&& element : elements(grid_view)) {
    if (gt != element.geometry().type()) {
      return false;
    }
  }
  return true;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH
