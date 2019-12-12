#ifndef DUNE_COPASI_DYNAMIC_POWER_LOCAL_BASIS_HH
#define DUNE_COPASI_DYNAMIC_POWER_LOCAL_BASIS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <vector>

namespace Dune::Copasi {

/**
 * @brief      This class describes dynamic power local basis.
 *
 * @tparam     Basis  The base local basis
 */
template<class Basis>
class DynamicPowerLocalBasis
{
public:
  using Traits = typename Basis::Traits;

  DynamicPowerLocalBasis(std::unique_ptr<Basis>&& basis, std::size_t power_size = 1)
    : _power_size(power_size)
    , _basis(basis)
  {
    assert(_power_size >= 0);
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  basis       The base local basis
   * @param[in]  power_size  The power size
   */
  DynamicPowerLocalBasis(const Basis& basis, std::size_t power_size = 1)
    : _power_size(power_size)
  {
    assert(_power_size >= 0);
    if constexpr (std::is_polymorphic_v<Basis>)
      _basis = std::unique_ptr<const Basis>(basis.clone());
    else
      _basis = std::make_unique<const Basis>(basis);
  }

  DynamicPowerLocalBasis(const DynamicPowerLocalBasis& other)
    : _power_size(other._power_size)
  {
    if constexpr (std::is_polymorphic_v<Basis>)
      _basis = std::unique_ptr<const Basis>(other._basis->clone());
    else
      _basis = std::make_unique<const Basis>(*other._basis);
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  power_size  The power size
   *
   * @tparam     <unnamed>   Internal use to allow default constructible base
   *                         local basis
   */
  template<
    bool default_constructible = std::is_default_constructible_v<Basis>,
    class = std::enable_if_t<default_constructible>>
  DynamicPowerLocalBasis(std::size_t power_size = 1)
    : DynamicPowerLocalBasis(std::make_unique<const Basis>(), power_size)
  {}

  /**
   * @brief      Returns the number of shape functions for this local finite
   *             element
   *
   * @return     The number of shape functions
   */
  unsigned int size() const { return _power_size * _basis->size(); }

  /**
   * @brief      Evaluation function for a given local coordinate
   *
   * @param[in]  in    The local coordinates
   * @param      out   The value of each shape function
   */
  inline void evaluateFunction(
    const typename Traits::DomainType& in,
    std::vector<typename Traits::RangeType>& out) const
  {
    auto f = [&](const auto& i, auto& o) { _basis->evaluateFunction(i, o); };
    populate_output(f, in, out);
  }

  /**
   * @brief      Jacobian evaluation function for a given local coordinate
   *
   * @param[in]  in    The local coordinates
   * @param      out   The value of the jacobian for each shape function
   */
  inline void evaluateJacobian(
    const typename Traits::DomainType& in,
    std::vector<typename Traits::JacobianType>& out) const
  {
    auto f = [&](const auto& i, auto& o) { _basis->evaluateJacobian(i, o); };
    populate_output(f, in, out);
  }

  /**
   * @brief      The partial derivate of the shape functions
   * @warning    This function is here only fulfill the interface requirements.
   *             This was never used nor tested.
   *
   * @param[in]  order  The order
   * @param[in]  in     The local coordinates
   * @param      out    The value of the partial derivate for each shape
   * function
   *
   * @tparam     dim    The grid entity dimension
   */
  void partial(const std::array<unsigned int, Traits::dimDomain>& order,
               const typename Traits::DomainType& in,
               std::vector<typename Traits::RangeType>& out) const
  {
    auto f = [&](const auto& i, auto& o) { _basis->partial(order, i, o); };
    populate_output(f, in, out);
  }

  /**
   * @brief      Polynomal order of the shape functions
   *
   * @return     The polynomal order of the shape functions
   */
  unsigned int order() const { return _basis->order(); }

private:
  /**
   * @brief      Helper class to polulate results for interface functions
   *
   * @param[in]  f     A function that evaluates a interface function for an
   *                   individual local basis function
   * @param[in]  in    The local coordinates
   * @param      out   The result of the evaluation of f over the power local
   *                   basis function
   *
   * @tparam     F     The function type
   * @tparam     In    The local coordinates type
   * @tparam     Out   The result populated type
   */
  template<class F, class In, class Out>
  inline void populate_output(const F& f, const In& in, Out& out) const
  {
    out.resize(_power_size * _basis->size());
    if (_power_size == 0)
      return;

    f(in, out);
    assert(out.size() == _basis->size());

    if (_power_size == 1)
      return;

    auto it = out.begin();
    auto copy_begin = it;

    // skip the first n values because they were already evaluated
    std::advance(it, _basis->size());
    auto copy_end = it;

    while (it != out.end()) {
      std::copy(copy_begin, copy_end, it);
      std::advance(it, _basis->size());
    }
  }

private:
  std::unique_ptr<const Basis> _basis;
  const std::size_t _power_size;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_POWER_LOCAL_BASIS_HH
