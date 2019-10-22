#ifndef DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH

#include <dune/localfunctions/common/localinterpolation.hh>

#include <dune/common/classname.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <iostream>

namespace Dune::Copasi {

/**
 * @brief      This class describes a dynamic power local interpolation.
 *
 * @tparam     Interpolation  The base interpolation
 */
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

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  interpolation  The base interpolation
   */
  DynamicPowerLocalInterpolation(const Interpolation& interpolation)
    : DynamicPowerLocalInterpolation(interpolation, 1)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  power_size  The power size
   *
   * @tparam     <unnamed>   Internal template argument to allow default
   *                         constructible base interpolations
   */
  template<
    class = std::enable_if_t<std::is_default_constructible_v<Interpolation>>>
  DynamicPowerLocalInterpolation(std::size_t power_size)
    : DynamicPowerLocalInterpolation(Interpolation{}, power_size)
  {}

  /**
   * @brief      Constructs a new instance.
   *
   * @tparam     <unnamed>   Internal template argument to allow default
   *                         constructible base interpolations
   */
  template<
    class = std::enable_if_t<std::is_default_constructible_v<Interpolation>>>
  DynamicPowerLocalInterpolation()
    : DynamicPowerLocalInterpolation(1)
  {}

  /**
   * @brief      Interpolation method
   *
   * @param[in]  f     Function to be interpolated. If its return type is a
   *                   dynamic vector, it has to have a size equal or greater
   *                   than the power size set to this ibject
   * @param      out   The vector of coefficients for the local finite element
   *
   * @tparam     F     The function type
   * @tparam     C     The coefficents type
   */
  template<typename F, typename C>
  void interpolate(const F& f, std::vector<C>& out) const
  {
    out.clear();

    static_assert(IsIndexable<typename F::RangeType>::value);
    constexpr bool dynamic_vector =
      IsIndexable<decltype(std::declval<typename F::RangeType>()[0])>::value;

    if (_power_size == 0)
      return;
    if constexpr (not dynamic_vector) {
      _interpolation.interpolate(f, out);
    } else {

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

        // resize out vector to the correct size
        if (i == 0)
          out.resize(base_out.size() * _power_size);

        // output iterator
        auto out_it = out.begin();

        // move output iterator to the current component to interpolate
        std::advance(out_it, i * base_out.size());

        // copy result into the output container
        std::copy(base_out.begin(), base_out.end(), out_it);
      }
    }
  }

private:
  Interpolation _interpolation;
  std::size_t _power_size;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
