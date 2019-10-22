#ifndef DUNE_COPASI_MULTIDOMAIN_ENTITY_TRANSFORMATION_HH
#define DUNE_COPASI_MULTIDOMAIN_ENTITY_TRANSFORMATION_HH

#include <dune/copasi/concepts/grid.hh>

namespace Dune::Copasi {

/**
 * @brief      Multidomain entity transformation
 * @details    This class transform an entity of a multidomain grid to a
 *             subdomain grid entity. Equivalent to the lambda function:
 * @code{.cpp}
 * auto entity_transform = [&](auto e){return
 * _grid->multiDomainEntity(entity);};
 * @endcode
 *
 * @tparam     Grid  The grid
 */
template<class Grid>
struct MultiDomainEntityTransformation
{
  static_assert(Concept::isMultiDomainGrid<Grid>());

  using MultiDomainEntity = typename Grid::Traits::template Codim<0>::Entity;
  using SubDomainEntity =
    typename Grid::SubDomainGrid::Traits::template Codim<0>::Entity;

  /**
   * @brief      Constructor.
   *
   * @param[in]  grid  A pointer to the grid
   */
  MultiDomainEntityTransformation(std::shared_ptr<const Grid> grid)
    : _grid(grid)
  {}

  /**
   * @brief      Function call operator.
   *
   * @param[in]  entity  The subdomain entity
   *
   * @return     The multidomain entity
   */
  const MultiDomainEntity& operator()(const SubDomainEntity& entity)
  {
    return _grid->multiDomainEntity(entity);
  }

private:
  const std::shared_ptr<const Grid> _grid;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MULTIDOMAIN_ENTITY_TRANSFORMATION_HH