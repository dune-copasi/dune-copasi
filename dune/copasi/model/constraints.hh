#ifndef DUNE_COPASI_MODEL_CONSTRAINTS_HH
#define DUNE_COPASI_MODEL_CONSTRAINTS_HH

#include <dune/copasi/grid/boundary_entity_mapper.hh>
#include <dune/copasi/model/functor_factory.hh>
#include <dune/copasi/model/local_domain.hh>

#include <dune/pdelab/basis/constraints/container_affine.hh>
#include <dune/pdelab/common/concurrency/shared_stash.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/common/parametertree.hh>

#include <function2/function2.hpp>

namespace Dune::Copasi {

/**
 * @brief Constraints parser and operator for a leaf basis
 * @warning Expressions with time dependency are not supported: the 'time'
 *          token will be evaluated to NaN
 *
 * @tparam GridView The grid view of the basis to constrain
 */
template<Dune::Concept::GridView GridView>
class Constraints : public TypeTree::LeafNode
{
  constexpr static auto dim = GridView::dimension;
  struct Data
  {
    std::unique_ptr<LocalDomain<dim>> local_domain;
    fu2::unique_function<FieldVector<double, 1>() const noexcept> constrain_fnc;
    Data(std::unique_ptr<LocalDomain<dim>> local_domain,
         fu2::unique_function<FieldVector<double, 1>() const noexcept> constrain_fnc)
      : local_domain{ std::move(local_domain) }
      , constrain_fnc{ std::move(constrain_fnc) }
    {
    }
  };

public:
  bool doConstrainBoundary() const { return _data_boundary and _data_boundary->constrain_fnc; }
  bool doConstrainSkeleton() const { return _data_skeleton and _data_skeleton->constrain_fnc; }
  bool doConstrainVolume() const { return _data_volume and _data_volume->constrain_fnc; }

  template<PDELab::Concept::MultiIndex MultiIndex, Dune::Concept::GridView EntitySet>
  using Container = PDELab::AffineConstraintsContainer<double, MultiIndex, EntitySet>;

  explicit Constraints(std::shared_ptr<BoundaryEntityMapper<GridView>> mapper,
                       const ParameterTree& config = {},
                       std::shared_ptr<const FunctorFactory<dim>> functor_factory = nullptr)
    : TypeTree::LeafNode()
    , _mapper{ mapper }
    , _data_volume([_config = config.sub("volume"),
                    _functor_factory = functor_factory]() -> std::unique_ptr<Data> {
      if (not _functor_factory)
        return nullptr;
      auto local_domain = std::make_unique<LocalDomain<dim>>();
      local_domain->in_volume = 1;
      // currently, time-dependent constraints are not possible in a generic way...
      local_domain->time = std::numeric_limits<double>::quiet_NaN();
      auto constraints =
        _functor_factory->make_scalar("constrain.volume", _config, *local_domain, 0);
      return std::make_unique<Data>(std::move(local_domain), std::move(constraints));
    })
    , _data_skeleton([_config = config.sub("skeleton"),
                      _functor_factory = functor_factory]() -> std::unique_ptr<Data> {
      if (not _functor_factory)
        return nullptr;
      auto local_domain = std::make_unique<LocalDomain<dim>>();
      // currently, time-dependent constraints are not possible in a generic way...
      local_domain->time = std::numeric_limits<double>::quiet_NaN();
      auto constraints =
        _functor_factory->make_scalar("constrain.skeleton", _config, *local_domain, 1);
      return std::make_unique<Data>(std::move(local_domain), std::move(constraints));
    })
    , _data_boundary([_config = config.sub("boundary"),
                      _functor_factory = functor_factory]() -> std::unique_ptr<Data> {
      if (not _functor_factory)
        return nullptr;
      auto local_domain = std::make_unique<LocalDomain<dim>>();
      // currently, time-dependent constraints are not possible in a generic way...
      local_domain->time = std::numeric_limits<double>::quiet_NaN();
      auto constraints =
        _functor_factory->make_scalar("constrain.boundary", _config, *local_domain, 1);
      return std::make_unique<Data>(std::move(local_domain), std::move(constraints));
    })
  {
  }

  void constrainVolume(const PDELab::Concept::LocalBasisLeaf auto& lbasis, auto& container)
  {
    if (lbasis.size() == 0)
      return;

    const auto& entity = lbasis.element();
    _data_volume->local_domain->entity_volume = entity.geometry().volume();

    const auto& lkeys = lbasis.finiteElement().localCoefficients();
    for (std::size_t dof = 0; dof != lbasis.size(); ++dof) {
      // the codim to which this dof is attached to
      unsigned int codim = lkeys.localKey(dof).codim();
      _data_volume->local_domain->position = entity.geometry().center();
      double constraint = _data_volume->constrain_fnc();
      if (constraint != std::numeric_limits<double>::max())
        container.addTranslationConstraint(lbasis.index(dof), constraint);
    }
  }

  void constrainBoundary(const Dune::Concept::Intersection auto& intersection,
                         const PDELab::Concept::LocalBasisLeaf auto& lbasis_in,
                         auto& container)
  {
    if (lbasis_in.size() == 0)
      return;

    add_intersection_constraints(1,
                                 _data_boundary,
                                 lbasis_in,
                                 intersection,
                                 intersection.geometryInInside(),
                                 intersection.indexInInside(),
                                 container);
  }

  void constrainSkeleton(const Dune::Concept::Intersection auto& intersection,
                         const PDELab::Concept::LocalBasisLeaf auto& lbasis_in,
                         const PDELab::Concept::LocalBasisLeaf auto& lbasis_out,
                         auto& container)
  {
    if (lbasis_in.size() == 0 and lbasis_out.size() == 0)
      return;

    if (lbasis_in.size() != 0)
      add_intersection_constraints(1,
                                   _data_skeleton,
                                   lbasis_in,
                                   intersection,
                                   intersection.geometryInInside(),
                                   intersection.indexInInside(),
                                   container);

    if (intersection.neighbor() and lbasis_out.size() != 0)
      add_intersection_constraints(-1,
                                   _data_skeleton,
                                   lbasis_out,
                                   intersection,
                                   intersection.geometryInOutside(),
                                   intersection.indexInOutside(),
                                   container);
  }

private:
  auto add_intersection_constraints(char dir_sign,
                                    const auto& data,
                                    const auto& lbasis,
                                    const auto& intersection,
                                    const auto& ig_geometry,
                                    auto face,
                                    auto& container) const
  {
    const auto& entity = lbasis.element();
    // find dof indices that belong to the intersection
    const auto& refelem = referenceElement(entity.geometry());
    data->local_domain->entity_volume = intersection.geometry().volume();
    data->local_domain->in_boundary = static_cast<double>(not intersection.neighbor());
    data->local_domain->in_skeleton = static_cast<double>(intersection.neighbor());

    const auto& lkeys = lbasis.finiteElement().localCoefficients();
    for (std::size_t dof = 0; dof != lbasis.size(); ++dof) {
      // the codim to which this dof is attached to
      unsigned int codim = lkeys.localKey(dof).codim();
      if (codim == 0)
        continue;
      // find the reference sub_entity index for this degree of freedom
      int sub_entity = lkeys.localKey(dof).subEntity();

      for (int j = 0; j != refelem.size(face, 1, codim); ++j) {
        if (sub_entity == refelem.subEntity(face, 1, j, codim)) {
          if (data->local_domain->in_boundary == _mapper->isBoundary(entity, sub_entity, codim)) {
            auto inside_pos = refelem.position(sub_entity, codim);
            data->local_domain->normal = dir_sign * intersection.unitOuterNormal(ig_geometry.local(inside_pos));
            data->local_domain->position = entity.geometry().global(inside_pos);
            double constraint = data->constrain_fnc();
            if (constraint != std::numeric_limits<double>::max())
              container.addTranslationConstraint(lbasis.index(dof), constraint);
          }
        }
      }
    }
  };

  std::shared_ptr<BoundaryEntityMapper<GridView>> _mapper;
  PDELab::SharedStash<Data> _data_volume, _data_skeleton, _data_boundary;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_CONSTRAINTS_HH
