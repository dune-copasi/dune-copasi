#ifndef DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_HH

#include <dune/copasi/dynamic_local_basis.hh>
#include <dune/copasi/dynamic_local_coefficients.hh>
#include <dune/copasi/dynamic_local_interpolation.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

template<class LocalFiniteElement>
class DynamicPowerLocalFiniteElement
{
  using BaseLocalBasis          = typename LocalFiniteElement::Traits::LocalBasisType;
  using BaseLocalCoefficients   = typename LocalFiniteElement::Traits::LocalCoefficientsType;
  using BaseLocalInterpolation  = typename LocalFiniteElement::Traits::LocalInterpolationType;

  using LocalBasis          = DynamicPowerLocalBasis<BaseLocalBasis>;
  using LocalCoefficients   = DynamicPowerLocalCoefficients<BaseLocalCoefficients>;
  using LocalInterpolation  = DynamicPowerLocalInterpolation<BaseLocalInterpolation>;

public:

  using Traits =  LocalFiniteElementTraits<LocalBasis,LocalCoefficients,LocalInterpolation>;

  DynamicPowerLocalFiniteElement
  ( const LocalFiniteElement& finite_element, 
    const LocalBasis& local_basis,
    const LocalCoefficients& local_coeficcients,
    const LocalInterpolation& local_interpolation
  ) : _finite_element(finite_element)
    , _basis(local_basis)
    , _coefficients(local_coeficcients)
    , _interpolation(local_interpolation)
  {}

  DynamicPowerLocalFiniteElement
  ( const LocalFiniteElement& finite_element, 
    const BaseLocalBasis& local_basis,
    const BaseLocalCoefficients& local_coeficcients,
    const BaseLocalInterpolation& local_interpolation,
    std::size_t power_size
  ) : _finite_element(finite_element)
    , _basis(local_basis,power_size)
    , _coefficients(local_coeficcients,power_size)
    , _interpolation(local_interpolation,power_size)
  {}

  DynamicPowerLocalFiniteElement (std::size_t power_size)
    : DynamicPowerLocalFiniteElement(LocalFiniteElement{},BaseLocalBasis{},BaseLocalCoefficients{},BaseLocalInterpolation{},power_size)
  {}

  const typename Traits::LocalBasisType& localBasis () const
  {
    return _basis;
  }

  const typename Traits::LocalCoefficientsType& localCoefficients () const
  {
    return _coefficients;
  }

  const typename Traits::LocalInterpolationType& localInterpolation () const
  {
    return _interpolation;
  }

  unsigned int size () const
  {
    return _basis.size();
  }

  GeometryType type () const
  {
    return _finite_element.type();
  }

private:

  LocalFiniteElement  _finite_element;
  LocalBasis          _basis;
  LocalCoefficients   _coefficients;
  LocalInterpolation  _interpolation;
};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_HH
