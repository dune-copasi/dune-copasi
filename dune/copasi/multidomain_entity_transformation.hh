#ifndef DUNE_COPASI_MULTIDOMAIN_ENTITY_TRANSFORMATION_HH
#define DUNE_COPASI_MULTIDOMAIN_ENTITY_TRANSFORMATION_HH

#include <dune/copasi/concepts/grid.hh>

namespace Dune::Copasi {

template<class Grid>
struct MultiDomainEntityTransformation
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using MultiDomainEntity = typename Grid::Traits::template Codim<0>::Entity;
  using SubDomainEntity =
    typename Grid::SubDomainGrid::Traits::template Codim<0>::Entity;

  MultiDomainEntityTransformation(std::shared_ptr<const Grid> grid)
    : _grid(grid)
  {}

  const MultiDomainEntity& operator()(const SubDomainEntity& entity)
  {
    return _grid->multiDomainEntity(entity);
  }

  const std::shared_ptr<const Grid> _grid;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MULTIDOMAIN_ENTITY_TRANSFORMATION_HH