#ifndef DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_HH

#include <dune/copasi/finite_element/dynamic_local_basis.hh>
#include <dune/copasi/finite_element/dynamic_local_coefficients.hh>
#include <dune/copasi/finite_element/dynamic_local_interpolation.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a dynamic power local finite element.
 *
 * @tparam     LocalFiniteElement  The base local finite element
 */
template<class LocalFiniteElement>
class DynamicPowerLocalFiniteElement
{
  using BaseLocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;
  using BaseLocalCoefficients =
    typename LocalFiniteElement::Traits::LocalCoefficientsType;
  using BaseLocalInterpolation =
    typename LocalFiniteElement::Traits::LocalInterpolationType;

  using LocalBasis = DynamicPowerLocalBasis<BaseLocalBasis>;
  using LocalCoefficients =
    DynamicPowerLocalCoefficients<BaseLocalCoefficients>;
  using LocalInterpolation =
    DynamicPowerLocalInterpolation<BaseLocalInterpolation>;

public:
  using Traits =
    LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  finite_element       The base local finite element
   * @param[in]  local_basis          The base local basis
   * @param[in]  local_coeficcients   The base local coeficcients
   * @param[in]  local_interpolation  The base local interpolation
   */
  DynamicPowerLocalFiniteElement(const LocalFiniteElement& finite_element,
                                 const LocalBasis& local_basis,
                                 const LocalCoefficients& local_coeficcients,
                                 const LocalInterpolation& local_interpolation)
    : _finite_element(finite_element)
    , _basis(local_basis)
    , _coefficients(local_coeficcients)
    , _interpolation(local_interpolation)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  finite_element  The base local finite element
   * @param[in]  power_size      The power size
   * @param[in]  local_basis          The base local basis
   * @param[in]  local_coeficcients   The base local coeficcients
   * @param[in]  local_interpolation  The base local interpolation
   */
  DynamicPowerLocalFiniteElement(const LocalFiniteElement& finite_element,
                                 std::size_t power_size)
    : _finite_element(finite_element)
    , _basis(_finite_element.localBasis(), power_size)
    , _coefficients(_finite_element.localCoefficients(), power_size)
    , _interpolation(_finite_element.localInterpolation(), power_size)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  power_size  The power size
   */
  DynamicPowerLocalFiniteElement(std::size_t power_size)
    : DynamicPowerLocalFiniteElement(LocalFiniteElement{}, power_size)
  {}

  /**
   * @brief      Returns the power local basis
   *
   * @return     The power local basis
   */
  const typename Traits::LocalBasisType& localBasis() const { return _basis; }

  /**
   * @brief      Returns the power local coefficents
   *
   * @return     The power local coefficients
   */
  const typename Traits::LocalCoefficientsType& localCoefficients() const
  {
    return _coefficients;
  }

  /**
   * @brief      Returns the power local interpolation
   *
   * @return     The power local interpolation
   */
  const typename Traits::LocalInterpolationType& localInterpolation() const
  {
    return _interpolation;
  }

  /**
   * @brief      Returns the number of degrees of freedom for this element
   *
   * @return     The size
   */
  unsigned int size() const { return _coefficients.size(); }

  /**
   * @brief      Returns the geometry type for this finite element
   *
   * @return     The geometry type
   */
  GeometryType type() const { return _finite_element.type(); }

private:
  LocalFiniteElement _finite_element;
  LocalBasis _basis;
  LocalCoefficients _coefficients;
  LocalInterpolation _interpolation;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_LOCAL_FINITE_ELEMENT_HH
