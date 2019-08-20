#ifndef DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH

#include <dune/localfunctions/common/localinterpolation.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

namespace Dune::Copasi {

template<class Interpolation>
class DynamicPowerLocalInterpolation
{
public:
  DynamicPowerLocalInterpolation(const Interpolation& interpolation,
                                 std::size_t power_size)
    : _interpolation(interpolation)
    , _power_size(power_size)
  {
    assert(_power_size >= 0);
  }

  DynamicPowerLocalInterpolation(const Interpolation& interpolation)
    : DynamicPowerLocalInterpolation(interpolation, 1)
  {}

  template<
    class = std::enable_if_t<std::is_default_constructible_v<Interpolation>>>
  DynamicPowerLocalInterpolation(std::size_t power_size)
    : DynamicPowerLocalInterpolation(Interpolation{}, power_size)
  {}

  template<
    class = std::enable_if_t<std::is_default_constructible_v<Interpolation>>>
  DynamicPowerLocalInterpolation()
    : DynamicPowerLocalInterpolation(1)
  {}

  template<typename F, typename C>
  void interpolate(const F& f, std::vector<C>& out) const
  {
    out.clear();

    if (_power_size == 0)
      return;

    if constexpr (std::is_same_v<typename F::RangeType,
                                 FieldVector<double, 1>>) {
      _interpolation.interpolate(f, out);
    } else if constexpr (std::is_same_v<typename F::RangeType,
                                        DynamicVector<double>>) {
      // output iterator
      auto out_it = out.begin();

      // convert f in callable
      auto&& callable =
        Impl::makeFunctionWithCallOperator<typename F::DomainType>(f);

      for (std::size_t i = 0; i < _power_size; ++i) {

        std::vector<C> base_out;

        // specializate callable for component i
        auto callable_i = [&](typename F::DomainType x) {
          return callable(x)[i];
        };

        // evaluate component i
        _interpolation.interpolate(callable_i, base_out);

        // copy result into the output container
        out.resize(base_out.size() * _power_size);
        std::copy(base_out.begin(), base_out.end(), out_it);

        // move output iterator to the next component to interpolate
        std::advance(out_it, base_out.size());
      }
    } else
      static_assert(AlwaysTrue<C>::value);
  }

private:
  Interpolation _interpolation;
  std::size_t _power_size;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
