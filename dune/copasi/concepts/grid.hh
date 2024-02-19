#ifndef DUNE_COPASI_CONCEPTS_GRID_HH
#define DUNE_COPASI_CONCEPTS_GRID_HH

#include <dune-copasi-config.hh>

#include <dune/grid/concepts/grid.hh>

/**
 * @ingroup Concepts
 */
namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

template<class E, class IS>
concept IndexableEntity =
  Dune::Concept::Entity<E> && Dune::Concept::IndexSet<IS> && requires(const E entity, const IS is) {
    {
      is.index(entity)
    } -> std::same_as<typename IS::IndexType>;
  };

/**
 * @brief   Concept for dune subdomain grids of multidomain grids
 * @details Checks whether the type fits the most of the dune interface
 *          for grid and is extended to a subdomain grid.
 *          Some checks are missing, but they are not important
 *          for the concept.
 */
template<class G>
concept SubDomainGrid =
  Dune::Concept::Grid<G> && std::same_as<G, typename G::MultiDomainGrid::SubDomainGrid> &&
  requires(const G cg, typename G::SubDomainIndex sub_domain) {
    {
      cg.multiDomainGrid()
    } -> std::convertible_to<const typename G::MultiDomainGrid&>;
    {
      cg.domain()
    } -> std::convertible_to<typename G::SubDomainIndex>;
    // Some other things missing but this is enough for now
  };

/**
 * @brief   Concept for dune multidomain grids
 * @details Checks whether the type fits the most of the dune interface
 *          for grid and is extended to a multidomain grid.
 *          Some checks are missing, but they are not important
 *          for the concept.
 */
template<class G>
concept MultiDomainGrid =
  SubDomainGrid<typename G::SubDomainGrid> && Dune::Concept::Grid<G> &&
  Dune::Concept::Grid<typename G::HostGrid> &&
  requires(const G cg, typename G::SubDomainIndex sub_domain) {
    {
      G::maxSubDomainIndexIsStatic()
    } -> std::convertible_to<bool>;
    typename G::LeafSubDomainInterfaceIterator;
    typename G::LevelSubDomainInterfaceIterator;
    typename G::LeafAllSubDomainInterfacesIterator;
    typename G::LevelAllSubDomainInterfacesIterator;
    {
      cg.maxSubDomainIndex()
    } -> std::convertible_to<typename G::SubDomainIndex>;
    {
      cg.subDomain(sub_domain)
    } -> std::convertible_to<const typename G::SubDomainGrid&>;
    {
      cg.maxAssignedSubDomainIndex()
    } -> std::convertible_to<typename G::SubDomainIndex>;
    {
      cg.supportLevelIndexSets()
    } -> std::convertible_to<bool>;
    requires requires(G g, const typename G::template Codim<0>::Entity& entity) {
      g.startSubDomainMarking();
      g.preUpdateSubDomains();
      g.updateSubDomains();
      g.postUpdateSubDomains();
      g.addToSubDomain(sub_domain, entity);
      g.removeFromSubDomain(sub_domain, entity);
      g.assignToSubDomain(sub_domain, entity);
      g.removeFromAllSubDomains(entity);
    };
  };

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_GRID_HH
