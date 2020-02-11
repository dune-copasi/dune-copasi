#ifndef DUNE_COPASI_GRID_MARK_STRIPES_HH
#define DUNE_COPASI_GRID_MARK_STRIPES_HH

#include <dune/grid/uggrid.hh>

#include <set>

namespace Dune::Copasi {

/**
 * @brief      Mark stripes for refinement
 * @details    If a cube element is surrounded by simplicies in oposite sides
 *             such cube is understood to be part of a stripe of cubes. In such
 *             case, this method should mark the grid such that a stripe is a stripe
 *             after refinment (e.g. refine in the direction of perpendicular to the simplices)
 * 
 * @warning    Grid type should be UGGrid
 * @warning    Many edge cases have not been tested, if you find a such a case, please report it as a bug
 * 
 * @param      grid         The grid
 * @param[in]  mark_others  If true, mark other type of entities for refinment.
 *
 * @tparam     Grid         The grid type
 */
template<class Grid>
void mark_stripes(Grid& grid, bool mark_others = true)
{
  constexpr std::size_t dim = Grid::dimension;
  using RuleType = typename Dune::UG_NS<dim>::RefinementRule;

  auto grid_view = grid.leafGridView();
  std::set<int> simplex_side;

  // Loop over the grid
  for (auto&& entity : Dune::elements(grid_view))
  {
    if (entity.type().isCube())
    {
      //register sides with simplices
      simplex_side.clear();
      for (auto&& ig : Dune::intersections(grid_view,entity))
        if(ig.neighbor())
          if (ig.outside().type().isSimplex())
            simplex_side.insert(ig.indexInInside()/2);

      // For now, lets only do it in 2D
      static_assert(dim == 2);

      if (simplex_side.size() == 1)
      {
        // side of the simplices
        int side = *(simplex_side.begin());
        // side to refine (permendicular)
        side = not side;
        // mark entity with a blue type refinment
        grid.mark(entity,RuleType::BLUE,side);
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
