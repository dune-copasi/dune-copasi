#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH

#include <dune/copasi/concepts/has_method.hh>
#include <dune/copasi/common/factory.hh>

#include <dune/localfunctions/lagrange/pk.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

template<class D, class R, int d, int k>
struct Factory<Dune::PkLocalFiniteElement<D,R,d,k>>
{
  template<class Context>
  static auto create(const Context& ctx)
  {
    if constexpr (Concept::has_method_order<Context>())
      assert(ctx.order() == k); // Order must be 0!
    static_assert(Concept::has_method_geometry_type<Context>());
    assert(ctx.geometry_type().isSimplex());
    return std::make_unique<Dune::PkLocalFiniteElement<D,R,d,k>>();
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH