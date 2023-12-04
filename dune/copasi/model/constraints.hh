#ifndef DUNE_COPASI_MODEL_CONSTRAINTS_HH
#define DUNE_COPASI_MODEL_CONSTRAINTS_HH

#include <dune/copasi/model/local_equations/functor_factory.hh>
#include <dune/copasi/model/local_equations/local_equations.hh>

#include <dune/pdelab/basis/constraints/container_affine.hh>
#include <dune/pdelab/common/concurrency/shared_stash.hh>
#include <dune/pdelab/concepts/multiindex.hh>

#include <dune/typetree/leafnode.hh>

#include <dune/common/parametertree.hh>

#include <function2/function2.hpp>

namespace Dune::Copasi {

template<Dune::Concept::Grid Grid>
class Constraints : public TypeTree::LeafNode
{
  static constexpr int dim = Grid::dimensionworld;
  using GridView = typename Grid::LeafGridView;

  struct Data
  {
    std::unique_ptr<LocalDomain<dim>> local_domain;
    fu2::unique_function<FieldVector<double, 1>() const noexcept> constrain_fnc;
  };

public:
  bool doConstrainBoundary() const { return _data_codim_1 and _data_codim_1->constrain_fnc; }
  bool doConstrainSkeleton() const { return _data_codim_1 and _data_codim_1->constrain_fnc; }
  bool doConstrainVolume() const { return _data_codim_0 and _data_codim_0->constrain_fnc; }

  template<PDELab::Concept::MultiIndex MultiIndex, Dune::Concept::GridView EntitySet>
  using Container = PDELab::AffineConstraintsContainer<double, MultiIndex, EntitySet>;

  explicit Constraints(const ParameterTree& config = {},
                       std::shared_ptr<const FunctorFactory<Grid>> functor_factory = nullptr)
    : TypeTree::LeafNode()
    , _data_codim_0(
        [_config = config, _functor_factory = functor_factory]() -> std::unique_ptr<Data> {
          if (not _functor_factory)
            return nullptr;
          auto local_domain = std::make_unique<LocalDomain<dim>>();
          local_domain->in_volume = 1;
          // currently, time-dependent constraints are not possible in a generic way...
          local_domain->time = std::numeric_limits<double>::quiet_NaN();
          auto constraints = _functor_factory->make_scalar("constrain", _config, *local_domain, 0);
          return std::make_unique<Data>(std::move(local_domain), std::move(constraints));
        })
    , _data_codim_1(
        [_config = config, _functor_factory = functor_factory]() -> std::unique_ptr<Data> {
          if (not _functor_factory)
            return nullptr;
          auto local_domain = std::make_unique<LocalDomain<dim>>();
          // currently, time-dependent constraints are not possible in a generic way...
          local_domain->time = std::numeric_limits<double>::quiet_NaN();
          auto constraints = _functor_factory->make_scalar("constrain", _config, *local_domain, 1);
          return std::make_unique<Data>(std::move(local_domain), std::move(constraints));
        })
  {
  }

  void constrainVolume(const PDELab::Concept::LocalBasisLeaf auto& lbasis, auto& container)
  {
    if (lbasis.size() == 0)
      return;

    const auto& entity = lbasis.element();
    _data_codim_0->local_domain->entity_volume = entity.geometry().volume();

    const auto& lkeys = lbasis.finiteElement().localCoefficients();
    for (std::size_t dof = 0; dof != lbasis.size(); ++dof) {
      // the codim to which this dof is attached to
      unsigned int codim = lkeys.localKey(dof).codim();
      if (codim != 0)
        continue;
      _data_codim_0->local_domain->position = entity.geometry().center();
      double constraint = _data_codim_0->constrain_fnc();
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

    const auto& entity_in = lbasis_in.element();
    // find dof indices that belong to the intersection
    const auto face = intersection.indexInInside();
    const auto& refelem = referenceElement(entity_in.geometry());
    _data_codim_1->local_domain->entity_volume = intersection.geometry().volume();
    _data_codim_1->local_domain->in_boundary = not intersection.neighbor();

    const auto& lkeys = lbasis_in.finiteElement().localCoefficients();
    for (std::size_t dof = 0; dof != lbasis_in.size(); ++dof) {
      // the codim to which this dof is attached to
      unsigned int codim = lkeys.localKey(dof).codim();
      if (codim == 0)
        continue;
      // find the reference sub_entity index for this degree of freedom
      int sub_entity = lkeys.localKey(dof).subEntity();
      for (int j = 0; j != refelem.size(face, 1, codim); ++j) {
        if (sub_entity == refelem.subEntity(face, 1, j, codim)) {
          auto inside_pos = refelem.position(sub_entity, codim);
          _data_codim_1->local_domain->normal =
            intersection.unitOuterNormal(intersection.geometryInInside().local(inside_pos));
          _data_codim_1->local_domain->position = entity_in.geometry().global(inside_pos);
          double constraint = _data_codim_1->constrain_fnc();
          if (constraint != std::numeric_limits<double>::max())
            container.addTranslationConstraint(lbasis_in.index(dof), constraint);
        }
      }
    }
    _data_codim_1->local_domain->in_boundary = false;
  }

  void constrainSkeleton(const Dune::Concept::Intersection auto& intersection,
                         const PDELab::Concept::LocalBasisLeaf auto& lbasis_in,
                         const PDELab::Concept::LocalBasisLeaf auto& lbasis_out,
                         auto& container)
  {
    if (lbasis_in.size() == 0 and lbasis_out.size() == 0)
      return;

    _data_codim_1->local_domain->in_skeleton = intersection.neighbor();
    constrainBoundary(intersection, lbasis_in, container);

    if (intersection.neighbor() and lbasis_out.size() != 0) {
      const auto& entity_out = lbasis_out.element();
      // find dof indices that belong to the intersection
      const auto face = intersection.indexInOutside();
      const auto& refelem = referenceElement(entity_out.geometry());
      _data_codim_1->local_domain->entity_volume = intersection.geometry().volume();

      const auto& lkeys = lbasis_out.finiteElement().localCoefficients();
      for (std::size_t dof = 0; dof != lbasis_out.size(); ++dof) {
        // the codim to which this dof is attached to
        unsigned int codim = lkeys.localKey(dof).codim();
        if (codim == 0)
          continue;
        // find the reference sub_entity index for this degree of freedom
        int sub_entity = lkeys.localKey(dof).subEntity();
        for (int j = 0; j != refelem.size(face, 1, codim); ++j) {
          if (sub_entity == refelem.subEntity(face, 1, j, codim)) {
            auto inside_pos = refelem.position(sub_entity, codim);
            _data_codim_1->local_domain->normal =
              intersection.unitOuterNormal(intersection.geometryInOutside().local(inside_pos));
            _data_codim_1->local_domain->position = entity_out.geometry().global(inside_pos);
            double constraint = _data_codim_1->constrain_fnc();
            if (constraint != std::numeric_limits<double>::max())
              container.addTranslationConstraint(lbasis_out.index(dof), constraint);
          }
        }
      }
    }
    _data_codim_1->local_domain->in_skeleton = false;
  }

private:
  PDELab::SharedStash<Data> _data_codim_0, _data_codim_1;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_CONSTRAINTS_HH
