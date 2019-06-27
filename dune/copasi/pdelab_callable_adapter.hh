#ifndef DUNE_COPASI_PDELAB_ADAPTERS_HH
#define DUNE_COPASI_PDELAB_ADAPTERS_HH

#include <utility>
#include <type_traits>

#include <dune/pdelab/common/function.hh>

#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

namespace Dune::Copasi {

template<class GV, class R>
class GridFunctionTraitsFromRange;

template<class GV, class RF, int dim>
class GridFunctionTraitsFromRange<GV,Dune::FieldVector<RF,dim>> 
  : public PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim>> {};

template<class GV, class RF>
class GridFunctionTraitsFromRange<GV,Dune::DynamicVector<RF>> 
  : public PDELab::GridFunctionTraits<GV,RF,-1,Dune::DynamicVector<RF>> {};


template<typename GV, typename R, typename F>
class GlobalCallableToGridFunctionAdapter
  : public Dune::PDELab::GridFunctionBase<GridFunctionTraitsFromRange<GV,R>,
                                          GlobalCallableToGridFunctionAdapter<GV,R,F> >
{
  GV gv;
  F f;
public:
  typedef GridFunctionTraitsFromRange<GV,R> Traits;

  //! construct from grid view
  GlobalCallableToGridFunctionAdapter (const GV& gv_, const F& f_) : gv(gv_), f(f_) {}

  //! get a reference to the grid view
  inline const GV& getGridView () const {return gv;}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xl,
                        typename Traits::RangeType& y) const
  {
    typename Traits::DomainType xg = e.geometry().global(xl);
    y = f(xg);
  }
};



template<typename GV, typename R, typename F>
class LocalCallableToGridFunctionAdapter
  : public Dune::PDELab::GridFunctionBase<GridFunctionTraitsFromRange<GV,R>,
                                          LocalCallableToGridFunctionAdapter<GV,R,F> >
{
  GV gv;
  F f;
public:
  typedef GridFunctionTraitsFromRange<GV,R> Traits;

  //! construct from grid view
  LocalCallableToGridFunctionAdapter (const GV& gv_, const F& f_) : gv(gv_), f(f_) {}

  //! get a reference to the grid view
  inline const GV& getGridView () const {return gv;}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xl,
                        typename Traits::RangeType& y) const
  {
    y = f(e,xl);
  }
};


/** \brief Create PDELab GridFunction from a callable f(x) that expects a global coordinate x */
template <typename GV, typename F>
auto makeGridFunctionFromCallable (const GV& gv, const F& f)
  -> typename std::enable_if<
    AlwaysTrue <
      decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>()))
      >::value,
    GlobalCallableToGridFunctionAdapter<
      GV,
        decltype(f(std::declval<typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate>())),
      F>
    >::type
{
  typedef typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate X;
  X x;
  typedef decltype(f(x)) Range;
  typedef GlobalCallableToGridFunctionAdapter<GV,Range,F> TheType;
  return TheType(gv,f);
}

/** \brief Create PDELab GridFunction from a callable f(e,x) that expects
    an entity e and a local coordinate x */
template <typename GV, typename F>
auto makeGridFunctionFromCallable (const GV& gv, const F& f)
  -> typename std::enable_if<
    AlwaysTrue <
      decltype(f(
                 std::declval<typename GV::template Codim<0>::Entity>(),
                 std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                 ))
      >::value,
    LocalCallableToGridFunctionAdapter<
      GV,
      decltype(f(
                 std::declval<typename GV::template Codim<0>::Entity>(),
                 std::declval<typename GV::template Codim<0>::Entity::Geometry::LocalCoordinate>()
                 )),
      F>
    >::type
{
  typedef typename GV::template Codim<0>::Entity E;
  E e;
  typedef typename E::Geometry::LocalCoordinate X;
  X x;
  typedef decltype(f(e,x)) Range;
  typedef LocalCallableToGridFunctionAdapter<GV,Range,F> TheType;
  return TheType(gv,f);
}


} // namespace Dune::Copasi

#endif // DUNE_COPASI_PDELAB_ADAPTERS_HH