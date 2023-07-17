#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_P0_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>

#include <dune/localfunctions/lagrange/p0.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

/**
 * @brief      Factory for P0LocalFiniteElement instances
 * @ingroup    Factory, FiniteElement
 * @tparam     <unnamed>  Template parameters of the P0LocalFiniteElement
 */
template<class D, class R, int d>
struct Factory<Dune::P0LocalFiniteElement<D,R,d>>
{
  /**
   * @brief      Create method
   *
   * @param      ctx   @ref DataContext containing the geometry type
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of P0LocalFiniteElement
   */
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