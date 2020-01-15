#ifndef DUNE_COPASI_GRID_MARK_STRIPES_HH
#define DUNE_COPASI_GRID_MARK_STRIPES_HH

#include <dune/grid/uggrid.hh>

#include <set>

namespace Dune::Copasi {



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
        // mark entity with a blute type refinment
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
