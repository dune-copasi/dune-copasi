#ifndef DUNE_COPASI_DYNAMIC_POWER_LOCAL_INTERPOLATION_HH
#define DUNE_COPASI_DYNAMIC_POWER_LOCAL_INTERPOLATION_HH

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

  DynamicPowerLocalInterpolation(std::unique_ptr<const Interpolation>&& interpolation,
                                 std::size_t power_size = 1)
    : _power_size(power_size)
    , _interpolation(std::move(interpolation))
  {
    assert(_power_size >= 0);
  }

  DynamicPowerLocalInterpolation(const Interpolation& interpolation,
                                 std::size_t power_size = 1)
    : _power_size(power_size)
  {
    assert(_power_size >= 0);
    if constexpr (std::is_polymorphic_v<Interpolation>)
      _interpolation = std::unique_ptr<const Interpolation>(interpolation.clone());
    else
      _interpolation = std::make_unique<const Interpolation>(interpolation);
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  power_size  The power size
   *
   * @tparam     <unnamed>   Internal template argument to allow default
   *                         constructible base interpolations
   */
  template<
    bool default_constructible = std::is_default_constructible_v<Interpolation>,
    class = std::enable_if_t<default_constructible>>
  DynamicPowerLocalInterpolation(std::size_t power_size = 1)
    : DynamicPowerLocalInterpolation(std::make_unique<Interpolation>(), power_size)
  {}

  /**
   * @brief      Copy constructor.
   *
   * @param[in]  other  The other interpolation to be copied
   *
   */
  DynamicPowerLocalInterpolation(const DynamicPowerLocalInterpolation& other)
    : _power_size(other._power_size)
  {
    if constexpr (std::is_polymorphic_v<Interpolation>)
      _interpolation = std::unique_ptr<const Interpolation>(other._interpolation->clone());
    else
      _interpolation = std::make_unique<const Interpolation>(*other._interpolation);
  }

  /**
   * @brief      Interpolation method
   *
   * @param[in]  f     Function to be interpolated. If its return type is a
   *                   dynamic vector, it has to have a size equal or greater
   *                   than the power size set to this object
   * @param      out   The vector of coefficients for the local finite element
   *
   * @tparam     F     The function type
   * @tparam     C     The coefficents type
   */
  template<typename F, typename C>
  void interpolate(const F& f, std::vector<C>& out) const
  {
    out.clear();

    // the range of the function must be indexable
    using Range = typename F::RangeType;
    static_assert(IsIndexable<Range>::value);

    // if field on the range is indexable, we assume they correspond to the 
    // components of a dynamic power finite element. 
    using RangeField = decltype(std::declval<Range>()[0]);
    
    constexpr bool dynamic_vector = IsIndexable<RangeField>::value;

    if (_power_size == 0)
    {
      out.resize(0);
      return;
    }

    if constexpr (not dynamic_vector) {
      assert(_power_size == 1);
      _interpolation->interpolate(f, out);
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
        _interpolation->interpolate(callable_i, base_out);

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
  const std::size_t _power_size;
  std::unique_ptr<const Interpolation> _interpolation;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_DYNAMIC_POWER_LOCAL_INTERPOLATION_HH
