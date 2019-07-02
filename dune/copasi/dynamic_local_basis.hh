#ifndef DUNE_COPASI_DYNAMIC_LOCAL_BASIS_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_BASIS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <iterator>
#include <vector>

namespace Dune::Copasi {

template<class Basis>
class DynamicPowerLocalBasis
{
public:
  using Traits = typename Basis::Traits;

  DynamicPowerLocalBasis(const Basis& basis, std::size_t power_size)
    : _basis(basis)
    , _size(power_size * _basis.size())
  {
    assert(power_size > 0);
  }

  DynamicPowerLocalBasis(const Basis& basis)
    : DynamicPowerLocalBasis(basis, 1)
  {}

  template<class = std::enable_if_t<std::is_default_constructible_v<Basis>>>
  DynamicPowerLocalBasis(std::size_t power_size)
    : DynamicPowerLocalBasis(Basis{}, power_size)
  {}

  template<class = std::enable_if_t<std::is_default_constructible_v<Basis>>>
  DynamicPowerLocalBasis()
    : DynamicPowerLocalBasis(1)
  {}

  unsigned int size() const { return _size; }

  inline void evaluateFunction(
    const typename Traits::DomainType& in,
    std::vector<typename Traits::RangeType>& out) const
  {
    auto f = [&](const auto& i, auto& o) { _basis.evaluateFunction(i, o); };
    populate_output(f, in, out);
  }

  inline void evaluateJacobian(
    const typename Traits::DomainType& in,
    std::vector<typename Traits::JacobianType>& out) const
  {
    auto f = [&](const auto& i, auto& o) { _basis.evaluateJacobian(i, o); };
    populate_output(f, in, out);
  }

  template<int dim>
  void partial(const std::array<unsigned int, dim>& order,
               const typename Traits::DomainType& in,
               std::vector<typename Traits::RangeType>& out) const
  {
    auto f = [&](const auto& i, auto& o) { _basis.partial(order, i, o); };
    populate_output(f, in, out);
  }

  unsigned int order() const { return _basis.order(); }

private:
  template<class F, class In, class Out>
  void populate_output(const F& f, const In& in, Out& out) const
  {
    out.clear();
    out.reserve(_size);
    f(in, out);

    assert(out.size() == _basis.size());
    out.resize(_size);

    auto it = out.begin();
    auto copy_begin = it;

    // skip the first n values because they were already evaluated
    std::advance(it, _basis.size());
    auto copy_end = it;

    while (it != out.end()) {
      std::copy(copy_begin, copy_end, it);
      std::advance(it, _basis.size());
    }
  }

private:
  Basis _basis;
  std::size_t _size;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_LOCAL_BASIS_HH
