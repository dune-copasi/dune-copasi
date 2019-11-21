#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH

#include <dune/copasi/concepts/has_method.hh>
#include <dune/copasi/common/factory.hh>

#include <dune/localfunctions/lagrange/p0.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

template<class D, class R, int d>
struct Factory<Dune::P0LocalFiniteElement<D,R,d>>
{
  template<class Context>
  static auto create(const Context& ctx)
  {
    if constexpr (Concept::has_method_order<Context>())
      assert(ctx.order() == 0); // Order must be 0!
    static_assert(Concept::has_method_geometry_type<Context>());
    return std::make_unique<Dune::P0LocalFiniteElement<D,R,d>>(ctx.geometry_type());
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH