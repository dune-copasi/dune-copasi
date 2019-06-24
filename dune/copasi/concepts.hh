#ifndef DUNE_COPASI_CONCEPTS_HH
#define DUNE_COPASI_CONCEPTS_HH

#include <dune/functions/common/functionconcepts.hh>

#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/common/function.hh>

#include <utility>

namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

struct PDELabGridFunction
{
  template<class GF>
  auto require(GF&& gf) -> decltype(
    requireBaseOf<Dune::PDELab::GridFunctionBase<typename GF::Traits,GF> >(gf)
  );
};


template<class GF>
static constexpr bool isPDELabGridFunction()
{ 
  return models<Concept::PDELabGridFunction, GF>(); 
}


template<class GV>
struct PDELabLocalCallable 
  : Refines<Dune::Functions::Concept::Callable<
              typename GV::template Codim<0>::Entity,
              typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate> 
            >
{
  using Entity = typename GV::template Codim<0>::Entity;
  using Domain = typename Entity::Geometry::LocalCoordinate;

  template<class F>
  auto require(F&& f) -> decltype(
    f(std::declval<Entity>(),std::declval<Domain>())
  );
};

template<class GV, class F>
static constexpr bool isPDELabLocalCallable()
{ 
  return models<Concept::PDELabLocalCallable<GV>, F>(); 
}

template<class GV>
struct PDELabGlobalCallable 
  : Refines<Dune::Functions::Concept::Callable<
              typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate> 
            >
{
  using Entity = typename GV::template Codim<0>::Entity;
  using Domain = typename Entity::Geometry::GlobalCoordinate;

  template<class F>
  auto require(F&& f) -> decltype(
    f(std::declval<Domain>())
  );
};

template<class GV, class F>
static constexpr bool isPDELabGlobalCallable()
{ 
  return models<Concept::PDELabGlobalCallable<GV>, F>(); 
}

template<class GV, class F>
static constexpr bool isPDELabCallable()
{ 
  return ( isPDELabLocalCallable<GV,F>() or isPDELabGlobalCallable<GV,F>() ); 
}


} // Dune::Copasi::Concept namespace

#endif // DUNE_COPASI_CONCEPTS_HH