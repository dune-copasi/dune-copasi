#ifndef DUNE_COPASI_CONCEPTS_GRID_HH
#define DUNE_COPASI_CONCEPTS_GRID_HH

#include <dune/functions/common/functionconcepts.hh>

#include <dune/geometry/type.hh>

/**
 * @ingroup Concepts
 */
namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

/**
 * @brief   Concept for dune grids
 * @details Checks whether the type fits the most of the dune interface 
 *          for grid. Some checks are missing, but they are not important 
 *          for the concept.
 */
struct Grid
{
  template<class G>
  auto require(G&& g) -> decltype(
    requireConvertible<int>(G::dimension),
    requireConvertible<int>(G::dimensionworld),
    requireType<typename G::LeafGridView>(),
    requireType<typename G::LevelGridView>(),
    // Codim missing!
    requireType<typename G::LeafIntersection>(),
    requireType<typename G::LevelIntersection>(),
    requireType<typename G::LeafIntersectionIterator>(),
    requireType<typename G::LevelIntersectionIterator>(),
    requireType<typename G::HierarchicIterator>(),
    requireType<typename G::LevelIndexSet>(),
    requireType<typename G::LeafIndexSet>(),
    requireType<typename G::GlobalIdSet>(),
    requireType<typename G::LocalIdSet>(),
    requireType<typename G::CollectiveCommunication>(),
    requireType<typename G::ctype>(),
    requireConvertible<int>(g.maxLevel()),
    requireConvertible<int>(g.size(std::declval<int>(), std::declval<int>())),
    requireConvertible<int>(g.size(std::declval<int>())),
    requireConvertible<int>(g.size(std::declval<int>(),
                                   std::declval<GeometryType>())),
    requireConvertible<int>(g.size(std::declval<GeometryType>())),
    requireConvertible<std::size_t>(g.numBoundarySegments()),
    requireConvertible<typename G::LevelGridView>(
      g.levelGridView(std::declval<int>())),
    requireConvertible<typename G::LeafGridView>(g.leafGridView()),
    requireConvertible<const typename G::GlobalIdSet&>(g.globalIdSet()),
    requireConvertible<const typename G::LocalIdSet&>(g.localIdSet()),
    requireConvertible<const typename G::LevelIndexSet&>(
      g.levelIndexSet(std::declval<int>())),
    requireConvertible<const typename G::LeafIndexSet&>(g.leafIndexSet()),
    g.globalRefine(std::declval<int>()),
    requireConvertible<bool>(
      g.mark(std::declval<int>(),
             std::declval<typename G::template Codim<0>::Entity>())),
    requireConvertible<int>(
      g.getMark(std::declval<typename G::template Codim<0>::Entity>())),
    requireConvertible<bool>(g.preAdapt()),
    requireConvertible<bool>(g.adapt()),
    g.postAdapt(),
    requireConvertible<typename G::CollectiveCommunication>(g.comm()),
    requireConvertible<std::size_t>(g.loadBalance())
    // entity missing!
  );
};

/**
 * @brief Check if a type is dune grid
 * 
 * @tparam G        The type to check
 * @return true     if the type is a dune grid
 * @return false    if the type is not a dune grid
 */
template<class G>
static constexpr bool
isGrid()
{
  return models<Concept::Grid, G>();
}

/**
 * @brief   Concept for dune multidomain grids
 * @details Checks whether the type fits the most of the dune interface 
 *          for grid and is extended to a multidomain grid. 
 *          Some checks are missing, but they are not important 
 *          for the concept.
 */
struct MultiDomainGrid : Refines<Dune::Copasi::Concept::Grid>
{
  template<class G>
  auto require(G&& g) -> decltype(
    requireType<typename G::HostGrid>(),
    requireType<typename G::SubDomainIndex>(),
    requireConvertible<typename G::SubDomainIndex>(g.maxSubDomainIndex()),
    requireConvertible<bool>(G::maxSubDomainIndexIsStatic()),
    requireType<typename G::SubDomainGrid>(),
    requireType<typename G::LeafSubDomainInterfaceIterator>(),
    requireType<typename G::LevelSubDomainInterfaceIterator>(),
    requireType<typename G::LeafAllSubDomainInterfacesIterator>(),
    requireType<typename G::LevelAllSubDomainInterfacesIterator>(),
    g.startSubDomainMarking(),
    g.preUpdateSubDomains(),
    g.updateSubDomains(),
    g.postUpdateSubDomains(),
    g.addToSubDomain(std::declval<typename G::SubDomainIndex>(),
                     std::declval<typename G::template Codim<0>::Entity>()),
    g.removeFromSubDomain(
      std::declval<typename G::SubDomainIndex>(),
      std::declval<typename G::template Codim<0>::Entity>()),
    g.assignToSubDomain(std::declval<typename G::SubDomainIndex>(),
                        std::declval<typename G::template Codim<0>::Entity>()),
    g.removeFromAllSubDomains(
      std::declval<typename G::template Codim<0>::Entity>()),
    requireConvertible<const typename G::SubDomainGrid&>(
      g.subDomain(std::declval<typename G::SubDomainIndex>())),
    requireConvertible<typename G::SubDomainIndex>(
      g.maxAssignedSubDomainIndex()),
    requireConvertible<bool>(g.supportLevelIndexSets()));
};


/**
 * @brief Check if a type is dune multidomain grid
 * 
 * @tparam G        The type to check
 * @return true     if the type is a dune multidomain grid
 * @return false    if the type is not a dune multidomain grid
 */
template<class G>
static constexpr bool
isMultiDomainGrid()
{
  return models<Concept::MultiDomainGrid, G>();
}

/**
 * @brief   Concept for dune subdomain grids of multidomain grids
 * @details Checks whether the type fits the most of the dune interface 
 *          for grid and is extended to a subdomain grid. 
 *          Some checks are missing, but they are not important 
 *          for the concept.
 */
struct SubDomainGrid : Refines<Dune::Copasi::Concept::Grid>
{
  template<class G>
  auto require(G&& g) -> decltype(
    requireType<typename G::HostGrid>(),
    requireType<typename G::MultiDomainGrid>(),
    requireConvertible<const typename G::MultiDomainGrid&>(g.multiDomainGrid()),
    requireConcept<Concept::MultiDomainGrid, typename G::MultiDomainGrid>(),
    requireConvertible<typename G::SubDomainIndex>(g.domain())
    // Some other things missing but this is enough
  );
};

/**
 * @brief Check if a type is dune subdomain grid
 * 
 * @tparam G        The type to check
 * @return true     if the type is a dune subdomain grid
 * @return false    if the type is not a dune subdomain grid
 */
template<class G>
static constexpr bool
isSubDomainGrid()
{
  return models<Concept::SubDomainGrid, G>();
}

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_GRID_HH