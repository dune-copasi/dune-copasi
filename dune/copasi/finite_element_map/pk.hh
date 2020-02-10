#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/common/data_context.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

namespace Dune::Copasi {

template<class GV, class DF, class RF, unsigned int k>
struct Factory<Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,k>>
{
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    using FEM = PDELab::PkLocalFiniteElementMap<GV,DF,RF,k>;
    static_assert(Ctx::has( Context::Tag<GV>{} ));
    return std::make_unique<FEM>(ctx.view( Context::Tag<GV>{} ));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH