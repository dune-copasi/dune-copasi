#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>
#include <dune/copasi/concepts/has_method.hh>
#include <dune/copasi/context/geometry_type.hh>


namespace Dune::Copasi {

template<class DF, class RF, int dim>
struct Factory<PDELab::P0LocalFiniteElementMap<DF,RF,dim>>
{
public:
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    using FEM = PDELab::P0LocalFiniteElementMap<DF,RF,dim>;
    if constexpr (Concept::has_method_geometry_type<Ctx>())
      return std::make_unique<FEM>(ctx.geometry_type());
    else if constexpr (Concept::has_method_grid_view<Ctx>())
    {
      if (not has_single_geometry_type(ctx.grid_view()))
        DUNE_THROW(InvalidStateException,"Grid view has to have only one grid view");

      GeometryType gt = ctx.grid_view().template begin<0>()->geometry().type();
      return std::make_unique<FEM>(gt);
    }
    else
      static_assert(AlwaysFalse<Ctx>::value, "Invalid provided context");
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH