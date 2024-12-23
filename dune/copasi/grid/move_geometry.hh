#ifndef DUNE_COPASI_GRID_MOVE_HH
#define DUNE_COPASI_GRID_MOVE_HH

#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/model/diffusion_reaction/local_equations.hh>

#include <dune/pdelab/common/local_container.hh>

#include <dune/grid/concepts/geometry.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/common/fvector.hh>

namespace Dune::Copasi {


namespace Impl {
template< class ct >
struct MultiLinearGeometrySimplexTraits : MultiLinearGeometryTraits<ct>
{
  template< int mydim, int cdim >
  struct CornerStorage
  {
    typedef std::array< FieldVector< ct, cdim > , mydim + 1> Type;
  };
};
}

/**
 * @brief Moves a geometry based on the local deformation equations
 * @pre geometry is a simplex
 */
template<Dune::Concept::Geometry Geometry>
Dune::Concept::Geometry auto
move_geometry(auto time,
              const PDELab::Concept::LocalBasis auto& lbasis,
              const PDELab::Concept::LocalConstContainer auto& lcoeff,
              auto& local_eqs,
              auto& fe_cache,
              const Geometry& geo,
              auto quad_proj)
{
  constexpr int mdim = Geometry::mydimension;
  constexpr int cdim = Geometry::coorddimension;

  assert(geo.type().isSimplex());
  assert(geo.corners() == mdim+1);

  local_eqs->time = time;
  local_eqs->in_volume = 1;
  local_eqs->entity_volume = geo.volume();

  std::array<typename Geometry::GlobalCoordinate, mdim+1> corner_storage;
  if (local_eqs->domain_deformation) {
    // get a quadrature rule with only corner points where we want to evaluate the deformation
    thread_local auto quad_rule = [](){
      QuadratureRule<typename Geometry::ctype, mdim> quad_rule;
      quad_rule.reserve(mdim+1);
      const auto& ref_el = ReferenceElements<typename Geometry::ctype, mdim>::general(GeometryTypes::simplex(mdim));
      for (int i = 0; i != ref_el.size(mdim); ++i)
        quad_rule.emplace_back(ref_el.position(i, mdim), std::numeric_limits<double>::quiet_NaN());
      return quad_rule;
    }();
    assert(quad_rule.size() == mdim+1);
    // evaluate deformation equations on each corner of the current geometry
    for (int i = 0; i != mdim+1; ++i) {
      // load trial functions into the equations cache
      forEachLeafNode(lbasis.tree(), [&](const auto& node) {
        if (node.size() == 0)
          return;
        auto& value = local_eqs->get_value(node);
        auto& gradient = local_eqs->get_gradient(node);
        value = 0.;
        // we have not created the geometry yet so we cannot evaluate the geometry jacobian
        // -> using the gradient would require a non-linear dependency on the geometry jacobian
        gradient = std::numeric_limits<std::decay_t<decltype(gradient)>>::quiet_NaN();
        fe_cache.bind(node.finiteElement(), quad_rule, quad_proj, true);
        const auto& phi = fe_cache.evaluateFunction(i);
        for (std::size_t dof = 0; dof != node.size(); ++dof)
          value += lcoeff(node, dof) * phi[dof];
      });

      // [[pre: numbering of reference element corners is the same as for the global geometry]]
      local_eqs->position = geo.corner(i);
      // evaluate deformation
      corner_storage[i] = local_eqs->position + local_eqs->domain_deformation();
    }
  } else {
    // copy corners of old geometry to new geometry
    for (int i = 0; i != mdim+1; ++i)
      corner_storage[i] = geo.corner(i);
  }

  using Traits = Impl::MultiLinearGeometrySimplexTraits<typename Geometry::ctype>;
  auto new_geo = MultiLinearGeometry<typename Geometry::ctype, mdim, cdim, Traits>(geo.type(), std::move(corner_storage));
  local_eqs->entity_volume = new_geo.volume();
  return new_geo;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MOVE_HH
