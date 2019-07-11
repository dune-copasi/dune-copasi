#ifndef DUNE_COPASI_CONCEPTS_HH
#define DUNE_COPASI_CONCEPTS_HH

#include <dune/functions/common/functionconcepts.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/function/callableadapter.hh>

#include <utility>

namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

struct PDELabGridFunction
{
  template<class GF>
  auto require(GF&& gf) -> decltype(
    requireBaseOf<Dune::PDELab::GridFunctionBase<typename GF::Traits, GF>>(gf));
};

template<class GF>
static constexpr bool
isPDELabGridFunction()
{
  return models<Concept::PDELabGridFunction, GF>();
}

template<class GV>
struct PDELabLocalCallable
  : Refines<Dune::Functions::Concept::Callable<
      typename GV::template Codim<0>::Entity,
      typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>>
{
  using Entity = typename GV::template Codim<0>::Entity;
  using Domain = typename Entity::Geometry::LocalCoordinate;

  template<class F>
  auto require(F&& f)
    -> decltype(f(std::declval<Entity>(), std::declval<Domain>()));
};

template<class GV, class F>
static constexpr bool
isPDELabLocalCallable()
{
  return models<Concept::PDELabLocalCallable<GV>, F>();
}

template<class GV>
struct PDELabGlobalCallable
  : Refines<Dune::Functions::Concept::Callable<
      typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>>
{
  using Entity = typename GV::template Codim<0>::Entity;
  using Domain = typename Entity::Geometry::GlobalCoordinate;

  template<class F>
  auto require(F&& f) -> decltype(f(std::declval<Domain>()));
};

template<class GV, class F>
static constexpr bool
isPDELabGlobalCallable()
{
  return models<Concept::PDELabGlobalCallable<GV>, F>();
}

template<class GV, class F>
static constexpr bool
isPDELabCallable()
{
  return (isPDELabLocalCallable<GV, F>() or isPDELabGlobalCallable<GV, F>());
}

struct Grid
{
  template<class G>
  auto require(G&& g) -> decltype(
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

template<class G>
static constexpr bool
isGrid()
{
  return models<Concept::Grid, G>();
}

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

template<class G>
static constexpr bool
isMultiDomainGrid()
{
  return models<Concept::MultiDomainGrid, G>();
}

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

template<class G>
static constexpr bool
isSubDomainGrid()
{
  return models<Concept::SubDomainGrid, G>();
}

struct TypeTree
{
  template<class T>
  auto require(T&& t) -> decltype(requireType<typename T::NodeTag>());
};

template<class T>
static constexpr bool
isTypeTree()
{
  return models<Concept::TypeTree, T>();
}

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_HH