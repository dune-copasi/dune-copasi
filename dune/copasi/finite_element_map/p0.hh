#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>
#include <dune/copasi/common/data_context.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

namespace Dune::Copasi {

template<class DF, class RF, int dim>
struct Factory<PDELab::P0LocalFiniteElementMap<DF,RF,dim>>
{
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    static_assert(Ctx::has( Context::Tag<GeometryType>{} ));

    using FEM = PDELab::P0LocalFiniteElementMap<DF,RF,dim>;
    return std::make_unique<FEM>(ctx.view( Context::Tag<GeometryType>{} ));
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH