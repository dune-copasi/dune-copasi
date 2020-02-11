#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>

#include <dune/localfunctions/lagrange/p0.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

template<class D, class R, int d>
struct Factory<Dune::P0LocalFiniteElement<D,R,d>>
{
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    using dCtx = std::decay_t<Ctx>;
    static_assert(dCtx::has( Context::Tag<GeometryType>{} ));
    const auto& gt = ctx.view( Context::Tag<GeometryType>{} );
    return std::make_unique<Dune::P0LocalFiniteElement<D,R,d>>(gt);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH