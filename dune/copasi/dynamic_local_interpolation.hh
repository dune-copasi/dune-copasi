#ifndef DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

namespace Dune::Copasi {

template<class F, class RF>
class LocalInterpolationWrapper
{
public:
  LocalInterpolationWrapper(const F& f, std::size_t power_size)
    : _f(f)
    , _power_size(power_size)
  {}

  void bind(std::size_t i) { _i = i; }

  template<class Domain, int dim>
  void evaluate(const Domain& x, FieldVector<RF, dim>& base_y) const
  {
    // evaluate f using a dynamic vector
    DynamicVector<RF> y(dim * _power_size);
    _f.evaluate(x, y);

    auto y_copy_begin = y.begin();
    std::advance(y_copy_begin, _i * dim);

    auto y_copy_end = y_copy_begin;
    std::advance(y_copy_end, dim);

    // copy sliced dynamic vector into a base output vector (e.g. field vector)
    std::copy(y_copy_begin, y_copy_end, base_y.begin());
  }

private:
  const F& _f;
  const std::size_t _power_size;
  std::size_t _i;
};

template<class Interpolation>
class DynamicPowerLocalInterpolation
{
public:
  DynamicPowerLocalInterpolation(const Interpolation& interpolation,
                                 std::size_t power_size)
    : _interpolation(interpolation)
    , _power_size(power_size)
  {
    assert(_power_size > 0);
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
    using RF = double;
    constexpr int range_dim = 1;

    // assume F is a dune-common function
    if constexpr (std::is_same_v<typename F::RangeType,
                                 FieldVector<RF, range_dim>>) {
      _interpolation.interpolate(f, out);
    } else if constexpr (std::is_same_v<typename F::RangeType,
                                        DynamicVector<RF>>) {
      // output iterator
      auto out_it = out.begin();

      // slice the output vector of f (dynamic) into an output vector
      // (field vector) that the interpolator actually expects
      LocalInterpolationWrapper<F, RF> wrapper(f, _power_size);

      for (int i = 0; i < _power_size; ++i) {
        std::vector<C> base_out;

        wrapper.bind(i);

        // evaluate component i
        _interpolation.interpolate(wrapper, base_out);

        // copy result into the output container
        out.reserve(base_out.size() * _power_size);
        std::copy(base_out.begin(), base_out.end(), out_it);

        // move output iterator to the next component to interpolate
        std::advance(out_it, base_out.size());
      }
    } else {
      static_assert(AlwaysTrue<F>::value, "Not known type");
    }
  }

private:
  std::size_t _power_size;
  Interpolation _interpolation;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
