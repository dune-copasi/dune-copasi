#ifndef DUNE_COPASI_GRID_BOUNDARY_ENTITY_MAPPER_HH
#define DUNE_COPASI_GRID_BOUNDARY_ENTITY_MAPPER_HH

#include <dune-copasi-config.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/grid/concepts/gridview.hh>

#include <memory>

namespace Dune::Copasi {


/**
 * @brief Mapper of boundary entities
 * @details Builds a map of entities with codim!=0 that are in the boundary.
 * This is helful to avoid referring to the same entities accessing sub-entities
 * of an intersection.
 * 
 * @tparam GridView The grid view to create the map from
 */
template<Dune::Concept::GridView GridView>
class BoundaryEntityMapper
{
public:
  BoundaryEntityMapper(GridView grid_view)
    : _mapper{ grid_view, [](GeometryType gt, int griddim) { return gt.dim() != static_cast<unsigned int>(griddim); } }
    , _boundary(_mapper.size(), false)
  {
    for (const auto& entity : elements(grid_view)) {
      for (const auto& intersection : intersections(grid_view, entity)) {
        if (intersection.boundary()) {
          auto face = intersection.indexInInside();
          const auto& refelem = referenceElement(entity.geometry());
          for (std::size_t codim = 1; codim != (GridView::dimension+1); ++codim) {
            unsigned int sz = refelem.size(face, 1, codim);
            for (unsigned int sub_entity = 0; sub_entity != sz; ++sub_entity) {
              unsigned int local_idx = refelem.subEntity(face, 1, sub_entity, codim);
              unsigned int idx = _mapper.subIndex(entity, local_idx, codim);
              _boundary[idx] = true;
            }
          }
        }
      }
    }
  }

  //! Whether a sub-entity is in the boundary of the grid view
  bool isBoundary(const auto& entity, auto local_index, auto codim) const
  {
    return _boundary[_mapper.subIndex(entity, local_index, codim)];
  }

  MultipleCodimMultipleGeomTypeMapper<GridView> _mapper;
  std::vector<bool> _boundary;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_BOUNDARY_ENTITY_MAPPER_HH
