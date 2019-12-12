#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH

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
    static_assert(Context::has(Signature::geometry_type));
    return std::make_unique<Dune::P0LocalFiniteElement<D,R,d>>(ctx.get(Signature::geometry_type));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH