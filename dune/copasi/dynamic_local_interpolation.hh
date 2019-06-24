#ifndef DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
#define DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

namespace Dune::Copasi {

template<class Interpolation>
class DynamicPowerLocalInterpolation
{
public:

  DynamicPowerLocalInterpolation(const Interpolation& interpolation, std::size_t power_size)
    : _interpolation(interpolation)
    , _power_size(power_size)
  {
    assert(_power_size > 0);
  }

  DynamicPowerLocalInterpolation(const Interpolation& interpolation)
    : DynamicPowerLocalInterpolation(interpolation,1)
  {}

  template<class = std::enable_if_t<std::is_default_constructible_v<Interpolation>>>
  DynamicPowerLocalInterpolation(std::size_t power_size)
    : DynamicPowerLocalInterpolation(Interpolation{},power_size)
  {}

  template<class = std::enable_if_t<std::is_default_constructible_v<Interpolation>>>
  DynamicPowerLocalInterpolation()
    : DynamicPowerLocalInterpolation(1)
  {}

  template<typename F, typename C>
  void interpolate (const F& f, std::vector<C>& out) const
  {
    using RF = double;
    constexpr int range_dim = 1;

    // assume F is a dune-common function
    if constexpr (std::is_same_v<typename F::RangeType,FieldVector<RF,range_dim>>)
    {
      _interpolation.interpolate(f,out);
    }
    else if constexpr (std::is_same_v<typename F::RangeType,DynamicVector<RF>>)
    {
      // output iterator
      auto out_it = out.begin();

      for (int i = 0; i < _power_size; ++i)
      {

        std::vector<C> base_out;

        // slice the output vector of f (dynamic) into an output vector 
        // (field vector) that the interpolator actually expects
        auto base_f = [&](auto x, auto base_y)
        {
          // evaluate f using a dynamic vector
          DynamicVector<RF> y(range_dim * _power_size);
          f.evaluate(x,y);
          
          auto y_copy_begin = y.begin();
          std::advance(y_copy_begin,i*range_dim);
          
          auto y_copy_end = y_copy_begin;
          std::advance(y_copy_end,range_dim+1);
          
          // copy sliced dynamic vector into a base output vector (e.g. field vector)
          std::copy(y_copy_begin,y_copy_end,base_y.begin());
        };

        // evaluate component i
        _interpolation.interpolate(base_f,base_out);

        // copy result into the output container
        out.reserve(base_out.size() * _power_size);
        std::copy(base_out.begin(),base_out.end(),out_it);

        // move output iterator to the next component to interpolate
        std::advance(out_it,base_out.size());
      }
    }
  }

private:

  std::size_t           _power_size;
  Interpolation         _interpolation;
};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_DYNAMIC_LOCAL_INTERPOLATION_HH
