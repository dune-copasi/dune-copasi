#ifndef DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH
#define DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/partitionset.hh>

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
template<class GridView, PartitionIteratorType partition = PartitionIteratorType::Interior_Partition>
std::enable_if_t<Capabilities::hasSingleGeometryType<typename GridView::Grid>::v,bool>
has_single_geometry_type(const GridView&)
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
template<class GridView, PartitionIteratorType partition = PartitionIteratorType::Interior_Partition>
std::enable_if_t<not Capabilities::hasSingleGeometryType<typename GridView::Grid>::v,bool>
has_single_geometry_type(const GridView& grid_view)
{
  auto it = grid_view.template begin<0,partition>();
  auto end = grid_view.template end<0,partition>();
  GeometryType gt;
  if (it != end)
    gt = it->geometry().type();
  for (; it != end; ++it)
    if (gt != it->geometry().type())
      return false;
  return true;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH
