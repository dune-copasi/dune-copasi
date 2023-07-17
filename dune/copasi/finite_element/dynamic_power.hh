#ifndef DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_HH
#define DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/finite_element/dynamic_power/local_basis.hh>
#include <dune/copasi/finite_element/dynamic_power/local_coefficients.hh>
#include <dune/copasi/finite_element/dynamic_power/local_interpolation.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

/**
 * @brief      This class describes a dynamic power local finite element.
 * @ingroup    FiniteElement
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
    DynamicPowerLocalInterpolation<BaseLocalInterpolation, typename LocalBasis::Traits::DomainType>;

public:
  using Traits =
    LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

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
    : _geometry_type(finite_element.type())
    , _basis(finite_element.localBasis(), power_size)
    , _coefficients(finite_element.localCoefficients(), power_size)
    , _interpolation(finite_element.localInterpolation(), power_size)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  power_size  The power size
   */
  template<
    bool default_constructible = std::is_default_constructible_v<LocalFiniteElement>,
    class = std::enable_if_t<default_constructible>>
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
   * @brief      Returns the power local coefficients
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
  GeometryType type() const { return _geometry_type; }

private:
  const GeometryType _geometry_type;
  const LocalBasis _basis;
  const LocalCoefficients _coefficients;
  const LocalInterpolation _interpolation;
};

/**
 * @brief      Factory for DynamicPowerLocalFiniteElement instances
 * @ingroup    Factory, FiniteElement
 * @tparam     BaseLocalFiniteElement  Base local finite element
 */
template<class BaseLocalFiniteElement>
struct Factory<DynamicPowerLocalFiniteElement<BaseLocalFiniteElement>>
{
public:
  /**
   * @brief      Create method
   *
   * @param      ctx   @ref DataContext containing the base finite element
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of DynamicPowerLocalFiniteElement
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    auto base_fe = Factory<BaseLocalFiniteElement>::create(std::forward<Ctx>(ctx));
    using FE = DynamicPowerLocalFiniteElement<BaseLocalFiniteElement>;
    return std::make_unique<FE>(*base_fe);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_POWER_LOCAL_FINITE_ELEMENT_HH
