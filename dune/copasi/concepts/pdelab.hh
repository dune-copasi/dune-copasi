#ifndef DUNE_COPASI_CONCEPTS_PDELAB_HH
#define DUNE_COPASI_CONCEPTS_PDELAB_HH

#include <dune/functions/common/functionconcepts.hh>

/**
 * @ingroup Concepts
 */
namespace Dune::Copasi::Concept {

using namespace Dune::Concept;

//! Concept for function traits in PDELab
struct PDELabFunctionTraits
{
  template<class F>
  auto require(F&& f)
    -> decltype(requireType<typename F::Traits>(),
                requireType<typename F::Traits::DomainFieldType>(),
                F::Traits::dimDomain,
                requireType<typename F::Traits::DomainType>(),
                requireType<typename F::Traits::RangeFieldType>(),
                F::Traits::dimRange,
                requireType<typename F::Traits::RangeType>());
};

//! Concept for functions in PDELab
struct PDELabFunction
{
  template<class F>
  auto require(F&& f) -> decltype(
    requireConcept<Dune::Copasi::Concept::PDELabFunctionTraits, F>(),
    f.evaluate(std::declval<const typename F::Traits::DomainType&>(),
               std::declval<typename F::Traits::RangeType&>()));
};

/**
 * @brief Check if a type is a function in PDELab
 *
 * @tparam GF       The type to check
 * @return true     if the type is a PDELab function
 * @return false    if the type is not a PDELab function
 */
template<class G>
static constexpr bool
isPDELabFunction()
{
  return models<Concept::PDELabFunction, G>();
}

//! Concept for grid function traits in PDELab
struct PDELabGridFunctionTraits : Refines<PDELabFunctionTraits>
{
  template<class GF>
  auto require(GF&& gf)
    -> decltype(requireType<typename GF::Traits::GridViewType>(),
                requireType<typename GF::Traits::ElementType>());
};

//! Concept for grid functions in PDELab
struct PDELabGridFunction
{
  template<class GF>
  auto require(GF&& gf) -> decltype(
    requireConcept<Dune::Copasi::Concept::PDELabGridFunctionTraits, GF>(),
    gf.evaluate(std::declval<const typename GF::Traits::DomainType&>(),
                std::declval<const typename GF::Traits::ElementType&>(),
                std::declval<typename GF::Traits::RangeType&>()),
    requireConvertible<typename GF::Traits::GridViewType>(gf.getGridView()));
};

/**
 * @brief Check if a type is a grid function in PDELab
 *
 * @tparam GF       The type to check
 * @return true     if the type is a PDELab grid function
 * @return false    if the type is not a PDELab grid function
 */
template<class GF>
static constexpr bool
isPDELabGridFunction()
{
  return models<Concept::PDELabGridFunction, GF>();
}

/**
 * @brief     Concept for local grid callables in PDELab
 * @warning   if the callable has two arguments and they do not fit the
 *            Entity and LocalCoordinate interface, this concept will fail
 *            inside the callable. Therefore, this concept an only be used
 *            when a local callable is actually required.
 *
 * @tparam GV   a grid view type
 */
template<class GV>
struct PDELabLocalCallable
{
  using Entity = typename GV::template Codim<0>::Entity;
  using Domain = typename Entity::Geometry::LocalCoordinate;

  template<class F>
  auto require(F&& f)
    -> decltype(f(std::declval<Entity>(), std::declval<Domain>()));
};

/**
 * @brief     Check if a type is a local callable usable as a function in PDELab
 * @warning   See restrictions on using this method on PDELabLocalCallable
 *
 * @tparam GF       The type to check
 * @tparam GV       The grid view type
 * @return true     if the type is a PDELab local callable
 * @return false    if the type is not a PDELab local callable
 */
template<class GV, class F>
static constexpr bool
isPDELabLocalCallable()
{
  return models<Concept::PDELabLocalCallable<GV>, F>();
}

/**
 * @brief     Concept for global grid callables in PDELab
 * @warning   if the callable has one argument and it does not fit the
 *            GlobalCoordinate interface, this concept will fail
 *            inside the callable. Therefore, this concept an only be used
 *            when a global callable is actually required.
 *
 * @tparam GV   a grid view type
 */
template<class GV>
struct PDELabGlobalCallable
{
  using Entity = typename GV::template Codim<0>::Entity;
  using Domain = typename Entity::Geometry::LocalCoordinate;

  template<class F>
  auto require(F&& f) -> decltype(f(std::declval<Domain>()));
};

/**
 * @brief     Check if a type is a global callable usable as a function in
 * PDELab
 * @warning   See restrictions on using this method on PDELabGlobalCallable
 *
 * @tparam GF       The type to check
 * @tparam GV       The grid view type
 * @return true     if the type is a PDELab global callable
 * @return false    if the type is not a PDELab global callable
 */
template<class GV, class F>
static constexpr bool
isPDELabGlobalCallable()
{
  return models<Concept::PDELabGlobalCallable<GV>, F>();
}

/**
 * @brief     Check if a type is a local or global callable usable as a function
 * in PDELab
 * @warning   See restrictions on using this method on PDELabLocalCallable and
 *            PDELabGlobalCallable
 *
 * @tparam GF       The type to check
 * @tparam GV       The type of a grid view
 * @return true     if the type is a PDELab callable
 * @return false    if the type is not a PDELab callable
 */
template<class GV, class F>
static constexpr bool
isPDELabCallable()
{
  return (isPDELabLocalCallable<GV, F>() or isPDELabGlobalCallable<GV, F>());
}

} // namespace Dune::Copasi::Concept

#endif // DUNE_COPASI_CONCEPTS_PDELAB_HH