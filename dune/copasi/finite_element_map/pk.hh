#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>
#include <dune/copasi/context/geometry_type.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

namespace Dune::Copasi {

template<class GV, class DF, class RF, int k>
struct Factory<PDELab::PkLocalFiniteElementMap<GV,DF,RF,k>>
{
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    using FEM = PDELab::PkLocalFiniteElementMap<GV,DF,RF,k>;
    static_assert(Ctx::has(Signature::grid_view));
    return std::make_unique<FEM>(ctx.get(Signature::grid_view));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_PK_HH