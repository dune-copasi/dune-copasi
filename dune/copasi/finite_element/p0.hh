#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/context/geometry_type.hh>

#include <dune/localfunctions/lagrange/p0.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

template<class D, class R, int d>
struct Factory<Dune::P0LocalFiniteElement<D,R,d>>
{
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    static_assert(Ctx::has( Context::Tag<GeometryType>{} ));
    const auto& gt = ctx.view( Context::Tag<GeometryType>{} );
    return std::make_unique<Dune::P0LocalFiniteElement<D,R,d>>(gt);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH