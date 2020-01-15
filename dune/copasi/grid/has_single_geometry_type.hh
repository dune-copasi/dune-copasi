#ifndef DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH
#define DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH

#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>

namespace Dune::Copasi {

template<class GridView>
std::enable_if_t<Capabilities::hasSingleGeometryType<typename GridView::Grid>::v,bool>
constexpr has_single_geometry_type(const GridView&)
{
  return true;
}

template<class GridView>
std::enable_if_t<not Capabilities::hasSingleGeometryType<typename GridView::Grid>::v,bool>
has_single_geometry_type(const GridView& grid_view)
{
  GeometryType gt = grid_view.template begin<0>()->geometry().type();
  for (auto&& element : elements(grid_view))
    if (gt != element.geometry().type())
    {
      std::cout << element.geometry().type() << std::endl;
      return false;
    }
  return true;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_HAS_SINGLE_GEOMETRY_TYPE_HH
