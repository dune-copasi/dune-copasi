#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

namespace Dune::Copasi {

/**
 * @brief      Factory for PkLocalFiniteElementMap instances
 * @ingroup    Factory, FiniteElementMap
 * @tparam     <unnamed>  Template parameters of the PkLocalFiniteElementMap
 */
template<class GV, class DF, class RF, unsigned int k>
struct Factory<Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,k>>
{
  /**
   * @brief      Create method
   *
   * @param      ctx   @ref DataContext containing a grid view of the type GV
   *
   * @tparam     Ctx   Universal reference to the @ref DataContext
   *
   * @return     Instance of PkLocalFiniteElementMap
   */
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    using dCtx = std::decay_t<Ctx>;
    using FEM = PDELab::PkLocalFiniteElementMap<GV,DF,RF,k>;
    static_assert(dCtx::has( Context::Tag<GV>{} ));
    return std::make_unique<FEM>(ctx.view( Context::Tag<GV>{} ));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH