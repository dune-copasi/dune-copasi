#ifndef DUNE_COPASI_GRID_MOVE_HH
#define DUNE_COPASI_GRID_MOVE_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/model/diffusion_reaction/local_equations.hh>

#include <dune/pdelab/common/local_container.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/geometry.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

namespace Dune::Copasi {

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
  if (not geo.affine())
    throw format_exception(NotImplemented{},
                           "Moving meshes with non-affine geometries are not yet implemented");

  local_eqs->time = time;
  local_eqs->in_volume = 1;

  constexpr int mdim = Geometry::mydimension;
  constexpr int cdim = Geometry::coorddimension;
  std::vector<typename Geometry::GlobalCoordinate> corner_storage(geo.corners());
  if (local_eqs->eval_position) {
    Dune::QuadraturePoint<typename Geometry::ctype, mdim> qp{
      FieldVector<typename Geometry::ctype, mdim>(0.), 0.
    };
    Dune::QuadratureRule<typename Geometry::ctype, Geometry::mydimension> quad_buffer;
    quad_buffer.resize(1, qp);
    for (int i = 0; i != geo.corners(); ++i) {
      quad_buffer[0] = { geo.local(geo.corner(i)), 0. };
      // evaluate concentrations at this corner
      forEachLeafNode(lbasis.tree(), [&](const auto& node) {
        if (node.size() == 0)
          return;
        auto& value = local_eqs->get_value(node);
        auto& gradient = local_eqs->get_gradient(node);
        value = 0.;
        // we have not created the geometry yet so we cannot evaluate the geometry jacobian
        // -> using the gradient would require a non-linear dependency on the geometry jacobian
        gradient = std::numeric_limits<std::decay_t<decltype(gradient)>>::quiet_NaN();
        fe_cache->bind(node.finiteElement(), quad_buffer, quad_proj, true);
        const auto& phi = fe_cache->evaluateFunction(0);
        for (std::size_t dof = 0; dof != node.size(); ++dof)
          value += lcoeff(node, dof) * phi[dof];
      });

      local_eqs->position = geo.corner(i);
      corner_storage[i] = local_eqs->eval_position();
    }
  } else {
    for (int i = 0; i != geo.corners(); ++i)
      corner_storage[i] = geo.corner(i);
  }
  // TODO change storage of std::vector to a reference!!!
  auto new_geo = MultiLinearGeometry<typename Geometry::ctype, mdim, cdim>(geo.type(), corner_storage);
  local_eqs->entity_volume = new_geo.volume();
  return new_geo;
}

template<class LBT, Dune::Concept::Grid Grid, class CellDataGridView>
void
move_grid(Grid& grid,
          const ParameterTree& config,
          auto time,
          const auto& basis,
          const auto& coefficients,
          std::shared_ptr<const FunctorFactory<Grid::dimensionworld>> functor_factory = nullptr,
          std::shared_ptr<const CellData<CellDataGridView, double>> grid_cell_data = nullptr)
{
  constexpr int dim = Grid::dimensionworld;
  MultipleCodimMultipleGeomTypeMapper mcmg{basis.entitySet(), mcmgVertexLayout()};
  using Coordinate = typename Grid::template Codim<dim>::Geometry::GlobalCoordinate;
  std::vector<Coordinate> positions(mcmg.size());
  auto lbasis = basis.localView();
  PDELab::LocalContainerBuffer lcontainer{ basis, &coefficients };
  auto fe_cache = std::make_unique<LocalBasisCache<LBT>>();
  auto local_eqs = DiffusionReaction::LocalEquations<dim>::make_mass(lbasis, config, functor_factory, grid_cell_data);
  // evaluate positions on each vertex
  for (const auto& entity : elements(basis.entitySet())) {
    lbasis.bind(entity);
    lcontainer.load(lbasis);
    const auto& geo = move_geometry(time, lbasis, lcontainer, local_eqs, fe_cache, entity.geometry(), std::identity{});
    for (int i = 0; i != entity.subEntities(dim); ++i) {
      // TODO are corners and sub-entity vertices equally enumerated?
      positions[mcmg.index(entity.template subEntity<dim>(i))] = geo.corner(i);
    }
  }
  // apply new vertex positions to grid
  for (const auto& vertex : vertices(basis.entitySet())) {
    auto pos = positions[mcmg.index(vertex)];
      if constexpr (requires {grid.setPosition(vertex, pos);} )
        grid.setPosition(vertex, pos);
      else if constexpr (requires {grid.hostGrid().setPosition(grid.hostEntity(vertex), pos);} )
        grid.hostGrid().setPosition(grid.hostEntity(vertex), pos);
      else
        throw format_exception(NotImplemented{}, "Moving of vertex positions is not implemented for this grid");
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MOVE_HH
