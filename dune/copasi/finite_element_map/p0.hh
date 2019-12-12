#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>

#include <dune/copasi/common/factory.hh>
#include <dune/copasi/grid/has_single_geometry_type.hh>
#include <dune/copasi/context/geometry_type.hh>


namespace Dune::Copasi {

template<class DF, class RF, int dim>
struct Factory<PDELab::P0LocalFiniteElementMap<DF,RF,dim>>
{
  template<class Ctx>
  static auto create(const Ctx& ctx)
  {
    using FEM = PDELab::P0LocalFiniteElementMap<DF,RF,dim>;
    if constexpr (Ctx::has(Signature::geometry_type))
      return std::make_unique<FEM>(ctx.get(Signature::geometry_type));
    else if constexpr (Ctx::has(Signature::grid_view))
    {
      const auto& gv = ctx.get(Signature::grid_view);
      std::cout << className<decltype(gv)>() << std::endl;
      if (not has_single_geometry_type(gv))
        DUNE_THROW(InvalidStateException,"Grid view has to have only one geometry type");

      GeometryType gt = gv.template begin<0>()->geometry().type();
      return std::make_unique<FEM>(gt);
    }
    else
      static_assert(AlwaysFalse<Ctx>::value, "Invalid Context");
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_MAP_P0_HH