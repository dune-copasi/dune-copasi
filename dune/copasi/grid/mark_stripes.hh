#ifndef DUNE_COPASI_GRID_MARK_STRIPES_HH
#define DUNE_COPASI_GRID_MARK_STRIPES_HH

#include <dune-copasi-config.h>

#include <dune/grid/uggrid.hh>

#include <list>

namespace Dune::Copasi {

/**
 * @brief      Mark stripes for refinement
 * @details    If a cube element is surrounded by non-cubes in opposite sides
 *             such cube is understood to be part of a stripe of cubes. In such
 *             case, this method should mark the grid such that a stripe is a
 *             stripe after refinement (e.g. refine in the direction of
 *             perpendicular to the non cubes)
 *
 * @param      grid         The grid
 * @param[in]  mark_others  If true, mark non-stripe entities for refinement.
 *
 * @tparam     dim         Dimension of the grid
 */
template<int dim>
void mark_stripes(UGGrid<dim>& grid, bool mark_others = true)
{
  using RuleType = typename UG_NS<dim>::RefinementRule;

  auto grid_view = grid.leafGridView();
  std::list<int> non_cube_side;

  // Loop over the grid
  for (auto&& entity : Dune::elements(grid_view))
  {
    if (entity.type().isCube())
    {
      // register side index with simplices (see DUNE cheatsheet for entity ids)
      non_cube_side.clear();
      for (auto&& ig : Dune::intersections(grid_view,entity))
        if(ig.neighbor())
          if (not ig.outside().type().isCube())
            non_cube_side.push_back(ig.indexInInside());

      bool is_stripe = false;

      // opposite facets have consecutive indexing (e.g. [2,3] are opposite)
      if (non_cube_side.size() == 2)
        is_stripe = !(non_cube_side.front()/2 - non_cube_side.back()/2);

      if (is_stripe)
      {
        // side orientation of the simplices
        [[maybe_unused]] int orientation = *(non_cube_side.begin())/2;

        if constexpr (dim == 2)
        {
          // mark entity with a blue type refinement
          grid.mark(entity,RuleType::BLUE,!(bool)orientation);
        }
        else if constexpr (dim == 3)
        {
          DUNE_THROW(NotImplemented,"\tStripes on 3D is not available yet!");
          // Need a mesh to correctly check which orientation needs which rule!
          // if (orientation == 0)
          //    grid.mark(entity,RuleType::HEX_BISECT_0_1);
          // if (orientation == 1)
          //    grid.mark(entity,RuleType::HEX_BISECT_0_2);
          // if (orientation == 2)
          //    grid.mark(entity,RuleType::HEX_BISECT_0_3);
        }
        else
        {
          DUNE_THROW(NotImplemented,
                     "\tStripe refinement not known for grids of dimension '"
                       << dim << "'");
        }
      }
      else if (mark_others)
      {
        grid.mark(1,entity);
      }
    }
    else if (mark_others)
    {
      grid.mark(1,entity);
    }
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MARK_STRIPES_HH
