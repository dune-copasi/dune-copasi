#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/model/local_equations/functor_factory.hh>
#include <dune/copasi/model/local_equations/local_equations.hh>

#include <dune/pdelab/common/concurrency/shared_stash.hh>
#include <dune/pdelab/common/execution.hh>
#include <dune/pdelab/common/tree_traversal.hh>
#include <dune/pdelab/operator/local_assembly/archetype.hh>
#include <dune/pdelab/operator/local_assembly/interface.hh>

#include <dune/grid/multidomaingrid/singlevalueset.hh>

#include <dune/geometry/type.hh>

#include <dune/common/fvector.hh>
#include <dune/common/overloadset.hh>

namespace Dune::Copasi {

enum class LocalOperatorType
{
  Stiffness,
  Mass
};

/**
 * @brief      This class describes a PDELab local operator for diffusion
 *             reaction systems.
 * @details    This class describre the operatrions for local integrals required
 *             for diffusion reaction system. The operator is only valid for
 *             entities contained in the entity set. The local finite element is
 *             used for caching shape function evaluations.
 *
 * @tparam     Basis    Basis
 * @tparam     LBT   Local basis traits
 */
template<PDELab::Concept::Basis TestBasis, class LBT, class ExecutionPolicy = PDELab::Execution::SequencedPolicy>
class LocalOperatorDiffusionReactionCG
{

  // utility types
  using DF = typename LBT::DomainFieldType;
  using RF = typename LBT::RangeFieldType;
  static constexpr int dim = LBT::dimDomain;

  using MembraneScalarFunction = typename LocalEquations<dim>::MembraneScalarFunction;
  using CompartmentNode = typename LocalEquations<dim>::CompartmentNode;
  struct Outflow
  {
    const MembraneScalarFunction& outflow;
    const CompartmentNode& source;
    Outflow(const MembraneScalarFunction& outflow, const CompartmentNode& source) : outflow{outflow}, source{source} {}
  };
  std::vector<Outflow> _outflow_i;
  std::vector<Outflow> _outflow_o;

  std::vector<std::size_t> _compartment2domain;

  TestBasis _test_basis;

  bool _is_linear = true;
  bool _has_outflow = true;

  PDELab::SharedStash<LocalBasisCache<LBT>> _fe_cache;
  PDELab::SharedStash<LocalEquations<dim>> _local_values_in;
  PDELab::SharedStash<LocalEquations<dim>> _local_values_out;

  ExecutionPolicy _execution_policy;

  template<class R, class X>
  struct PseudoJacobian
  {
    void accumulate(const auto& ltest,
                    auto test_dof,
                    const auto& ltrial,
                    auto trial_dof,
                    auto value)
    {
      _r.accumulate(ltest, test_dof, _z(ltrial, trial_dof) * value);
    }

    R& _r;
    const X& _z;
  };

  template<class R, class X>
  PseudoJacobian(R&, const X&) -> PseudoJacobian<R, X>;

  // mono-domains contain all entities, thus all domain sets are the same
  struct MonoDomainSet
  {
    static constexpr std::true_type contains(auto) { return {}; }
    friend constexpr std::true_type operator==(const MonoDomainSet&, const MonoDomainSet&)
    {
      return {};
    }
  };

  auto subDomains(const Dune::Concept::Entity auto& entity) const
  {
    if constexpr (Concept::MultiDomainGrid<typename TestBasis::EntitySet::Grid>)
      return _test_basis.entitySet().indexSet().subDomains(entity);
    else
      return MonoDomainSet{};
  }

  std::size_t integrationOrder(const auto&... lbasis_pack) const {
    std::size_t order = 0;
    std::size_t intorderadd = 0;
    std::size_t quad_factor = 2;
    Hybrid::forEach(std::tie(lbasis_pack...), [&](const auto& lbasis){
      forEachLeafNode(lbasis.tree(), [&](const auto& node) {
        if (node.size() != 0)
          order = std::max<std::size_t>(order, node.finiteElement().localBasis().order());
      });
    });
    return intorderadd + quad_factor * order;
  }

public:
  constexpr static std::true_type localAssembleDoVolume() noexcept { return {}; }

  constexpr static auto localAssembleDoSkeleton() noexcept
  {
    return std::bool_constant<Concept::MultiDomainGrid<typename TestBasis::EntitySet::Grid>>{};
  }

  auto executionPolicy() const {
    return _execution_policy;
  }

  bool localAssembleDoBoundary() const noexcept { return _has_outflow; }

  bool localAssembleSkipIntersection(const Dune::Concept::Intersection auto& intersection) const noexcept
  {
    // boundary case
    if (not intersection.neighbor())
      return false;

    // if domain sets are the same this is not an domain interface and should be
    // skipped
    return (subDomains(intersection.inside()) == subDomains(intersection.outside()));
  }

  bool localAssembleIsLinear() const noexcept { return _is_linear; }

  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  basis
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   */
  LocalOperatorDiffusionReactionCG(const PDELab::Concept::Basis auto& test_basis,
                                   LocalOperatorType lop_type,
                                   bool is_linear,
                                   const ParameterTree& config,
                                   std::shared_ptr<const FunctorFactory<dim>> functor_factory,
                                   ExecutionPolicy execution_policy = {})
    : _test_basis{ test_basis }
    , _is_linear{ is_linear }
    , _fe_cache([]() { return std::make_unique<LocalBasisCache<LBT>>(); })
    , _local_values_in([_lop_type = lop_type,
                     _basis = _test_basis,
                     _config = config,
                     _functor_factory = functor_factory]() {
      std::unique_ptr<LocalEquations<dim>> ptr;
      if (_lop_type == LocalOperatorType::Mass)
        ptr = LocalEquations<dim>::make_mass(_basis.localView(), _config, _functor_factory);
      else if (_lop_type == LocalOperatorType::Stiffness)
        ptr = LocalEquations<dim>::make_stiffness(_basis.localView(), _config, _functor_factory);
      if (not ptr)
        std::terminate();
      return ptr;
    })
    , _local_values_out([_lop_type = lop_type,
                         _basis = _test_basis,
                         _config = config,
                         _functor_factory = std::move(functor_factory)]() {
      std::unique_ptr<LocalEquations<dim>> ptr;
      if (_lop_type == LocalOperatorType::Mass)
        ptr = LocalEquations<dim>::make_mass(_basis.localView(), _config, _functor_factory);
      else if (_lop_type == LocalOperatorType::Stiffness)
        ptr = LocalEquations<dim>::make_stiffness(_basis.localView(), _config, _functor_factory);
      if (not ptr)
        std::terminate();
      return ptr;
    })
    , _execution_policy{ execution_policy }
  {
    auto lbasis = _test_basis.localView();
    if (_test_basis.entitySet().size(0) == 0)
      return;
    lbasis.bind(*_test_basis.entitySet().template begin<0>());
    _has_outflow = false;
    forEachLeafNode(lbasis.tree(), [&](const auto& ltrial_node) {
      const auto& eq = _local_values_in->get_equation(ltrial_node);
      _has_outflow |= not eq.outflow.empty();
    });
    lbasis.unbind();

    if constexpr (Concept::MultiDomainGrid<typename TestBasis::EntitySet::Grid>)
      forEachNode(lbasis.tree(),
                  overload(
                    [&](const Concept::CompartmentLocalBasisNode auto& /*ltrial_node*/, auto path) {
                      auto compartment = back(path);
                      _compartment2domain.resize(compartment + 1);
                      _compartment2domain[compartment] =
                        _test_basis.subSpace(path).entitySet().grid().domain();
                    },
                    [&](const auto& /*ltrial_node*/) {}));
    else
      _compartment2domain.assign(1, std::numeric_limits<std::size_t>::max());
  }

  /**
   * @brief      Pattern volume
   * @details    This method links degrees of freedom between trial and test
   *             basiss taking into account the structure of the reaction term
   *
   * @param[in]  lfsu          The trial local function basis
   * @param[in]  lfsv          The test local function basis
   * @param      pattern       The local pattern
   *
   * @tparam     LFSU          The trial local function basis
   * @tparam     LFSV          The test local function basis
   * @tparam     LocalPattern  The local pattern
   */
  void localAssemblePatternVolume(const PDELab::Concept::LocalBasis auto& ltrial,
                      const PDELab::Concept::LocalBasis auto& ltest,
                      auto& lpattern) const noexcept
  {
    forEachLeafNode(ltest.tree(), [&](const auto& ltest_node) {
      const auto& ltrial_node = PDELab::containerEntry(ltrial.tree(), ltest_node.path());
      const auto& eq = _local_values_in->get_equation(ltrial_node);
      if (eq.reaction) {
        for (const auto& jacobian_entry : eq.reaction.compartment_jacobian) {
          const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
              lpattern.addLink(ltest_node, dof_i, wrt_lbasis, dof_j);
        }
      }

      if (eq.storage) {
        for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
          for (std::size_t dof_j = 0; dof_j != ltrial_node.size(); ++dof_j)
            lpattern.addLink(ltest_node, dof_i, ltrial_node, dof_j);
        for (const auto& jacobian_entry : eq.storage.compartment_jacobian) {
          const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
              lpattern.addLink(ltest_node, dof_i, wrt_lbasis, dof_j);
        }
      }

      if (eq.velocity) {
        for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
          for (std::size_t dof_j = 0; dof_j != ltrial_node.size(); ++dof_j)
            lpattern.addLink(ltest_node, dof_i, ltrial_node, dof_j);
        for (const auto& jacobian_entry : eq.velocity.compartment_jacobian) {
          const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
              lpattern.addLink(ltest_node, dof_i, wrt_lbasis, dof_j);
        }
      }

      for (const auto& diffusion : eq.cross_diffusion) {
        const auto& wrt_lbasis = diffusion.wrt.to_local_basis_node(ltrial);

        for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
          for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
            lpattern.addLink(ltest_node, dof_i, wrt_lbasis, dof_j);

        for (const auto& jacobian_entry : diffusion.compartment_jacobian) {
          const auto& jac_wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j != jac_wrt_lbasis.size(); ++dof_j)
              lpattern.addLink(ltest_node, dof_i, jac_wrt_lbasis, dof_j);
        }
      }
    });
  }

  void localAssemblePatternSkeleton(const Dune::Concept::Intersection auto& intersection,
                                    const PDELab::Concept::LocalBasis auto& ltrial_in,
                                    const PDELab::Concept::LocalBasis auto& ltest_in,
                                    const PDELab::Concept::LocalBasis auto& ltrial_out,
                                    const PDELab::Concept::LocalBasis auto& ltest_out,
                                    auto& lpattern_in_in,
                                    auto& lpattern_in_out,
                                    auto& lpattern_out_in,
                                    auto& lpattern_out_out) noexcept
  {
    forEachLeafNode(ltest_in.tree(), [&](const auto& ltest_node_in) {
      if (ltest_node_in.size() == 0)
        return;
      const auto& ltrial_node = PDELab::containerEntry(ltrial_in.tree(), ltest_node_in.path());
      const auto& eq = _local_values_in->get_equation(ltrial_node);
      for (const auto& outflow_i : eq.outflow) {
        for (const auto& jacobian_entry : outflow_i.compartment_jacobian) {
          const auto& wrt_lbasis_in = jacobian_entry.wrt.to_local_basis_node(ltrial_in);
          const auto& wrt_lbasis_out = jacobian_entry.wrt.to_local_basis_node(ltrial_out);
          for (std::size_t dof_i = 0; dof_i != ltest_node_in.size(); ++dof_i) {
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis_in.size(); ++dof_j)
              lpattern_in_in.addLink(ltest_node_in, dof_i, wrt_lbasis_in, dof_j);
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis_out.size(); ++dof_j)
              lpattern_in_out.addLink(ltest_node_in, dof_i, wrt_lbasis_out, dof_j);
          }
        }
      }
    });

    if (not intersection.neighbor())
      return;

    forEachLeafNode(ltest_out.tree(), [&](const auto& ltest_node_out) {
      if (ltest_node_out.size() == 0)
        return;
      const auto& ltrial_node = PDELab::containerEntry(ltrial_out.tree(), ltest_node_out.path());
      const auto& eq = _local_values_in->get_equation(ltrial_node);
      for (const auto& outflow_o : eq.outflow) {
        for (const auto& jacobian_entry : outflow_o.compartment_jacobian) {
          const auto& wrt_lbasis_in = jacobian_entry.wrt.to_local_basis_node(ltrial_in);
          const auto& wrt_lbasis_out = jacobian_entry.wrt.to_local_basis_node(ltrial_out);
          for (std::size_t dof_i = 0; dof_i != ltest_node_out.size(); ++dof_i) {
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis_in.size(); ++dof_j)
              lpattern_out_in.addLink(ltest_node_out, dof_i, wrt_lbasis_in, dof_j);
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis_out.size(); ++dof_j)
              lpattern_out_out.addLink(ltest_node_out, dof_i, wrt_lbasis_out, dof_j);
          }
        }
      }
    });

  }

  void localAssemblePatternBoundary(const Dune::Concept::Intersection auto& intersection,
                                    const PDELab::Concept::LocalBasis auto& ltrial_in,
                                    const PDELab::Concept::LocalBasis auto& ltest_in,
                                    auto& lpattern_in) noexcept
  {
    localAssemblePatternSkeleton(intersection, ltrial_in, ltest_in, ltrial_in, ltest_in, lpattern_in, lpattern_in, lpattern_in, lpattern_in);
  }

  /**
   * @brief      The volume integral
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function basis
   * @param[in]  x     The local coefficient vector
   * @param[in]  lfsv  The test local function basis
   * @param      r     The local residual vector
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function basis
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function basis
   * @tparam     R     The local residual vector
   */

  void localAssembleVolume(auto time,
                           const PDELab::Concept::LocalBasis auto& ltrial,
                           const PDELab::Concept::LocalConstContainer auto& lcoefficients,
                           const PDELab::Concept::LocalBasis auto& ltest,
                           PDELab::Concept::LocalMutableContainer auto& lresidual) noexcept
  {
    if (ltrial.size() == 0)
      return;

    const auto& entity = ltrial.element();
    const auto& geo = entity.geometry();
    using Geometry = std::decay_t<decltype(geo)>;
    std::optional<typename Geometry::JacobianInverse> geojacinv_opt;

    _local_values_in->time = time;
    _local_values_in->entity_volume = geo.volume();
    _local_values_in->in_volume = 1;

    // Call the gridContext to update the gmsh_id
    _local_values_in->update_gmsh_id(entity);
    // Call the gridContext to update the grid variables
    _local_values_in->update_grid_data(entity);

    auto intorder = integrationOrder(ltrial);
    const auto& quad_rule = QuadratureRules<DF, dim>::rule(geo.type(), intorder);

    // loop over quadrature points
    for (std::size_t q = 0; q != quad_rule.size(); ++q) {
      const auto [position, weight] = quad_rule[q];
      if (not geojacinv_opt or not geo.affine())
        geojacinv_opt.emplace(geo.jacobianInverse(position));
      const auto& geojacinv = *geojacinv_opt;
      _local_values_in->position = geo.global(position);
      auto factor = weight * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      forEachLeafNode(ltrial.tree(), [&](const auto& node) {
        if (node.size() == 0)
          return;
        auto& value = _local_values_in->get_value(node);
        auto& gradient = _local_values_in->get_gradient(node);
        value = 0.;
        gradient = 0.;
        _fe_cache->bind(node.finiteElement(), quad_rule);
        const auto& phi = _fe_cache->evaluateFunction(q);
        const auto& jacphi = _fe_cache->evaluateJacobian(q);
        for (std::size_t dof = 0; dof != node.size(); ++dof) {
            value     += lcoefficients(node, dof) * phi[dof];
            gradient  += lcoefficients(node, dof) * (jacphi[dof] * geojacinv)[0];
        }
      });

      // contribution for each component
      forEachLeafNode(ltest.tree(), [&](const auto& ltest_node) {
        if (ltest_node.size() == 0)
          return;
        const auto& eq =
          _local_values_in->get_equation(PDELab::containerEntry(ltrial.tree(), ltest_node.path()));

        _fe_cache->bind(ltest_node.finiteElement(), quad_rule);
        const auto& psi = _fe_cache->evaluateFunction(q);
        const auto& jacpsi = _fe_cache->evaluateJacobian(q);

        RF scalar = eq.reaction ? RF{-eq.reaction()} : 0.;
        scalar += eq.storage ? RF{eq.value * eq.storage()} : 0.;
        auto flux = eq.velocity ? eq.velocity() * eq.value[0] : FieldVector<RF,dim>(0.);
        for (const auto& diffusion : eq.cross_diffusion)
          flux -= diffusion(diffusion.wrt.gradient);
        for (std::size_t dof = 0; dof != ltest_node.size(); ++dof)
          lresidual.accumulate(ltest_node, dof, (scalar * psi[dof] -dot(flux, (jacpsi[dof] * geojacinv)[0] )) * factor);
      });
    }

    _local_values_in->clear();
    _local_values_in->in_volume = 0;
  }

  /**
   * @brief      The jacobian volume integral for matrix free operations
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function basis
   * @param[in]  x     The local coefficient vector
   * @param[in]  z     The local position in the trial basis to which to apply
   *                   the Jacobian.
   * @param[in]  lfsv  The test local function basis
   * @param      r     The resulting vector
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function basis
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function basis
   * @tparam     R     The resulting vector
   */
  void localAssembleJacobianVolumeApply(
    auto time,
    const PDELab::Concept::LocalBasis auto& ltrial,
    const PDELab::Concept::LocalConstContainer auto& llin_point,
    const PDELab::Concept::LocalConstContainer auto& lapp_point,
    const PDELab::Concept::LocalBasis auto& ltest,
    PDELab::Concept::LocalMutableContainer auto& ljacobian) noexcept

  {
    PseudoJacobian mat{ ljacobian, lapp_point };
    if (localAssembleIsLinear())
      localAssembleVolume(time, ltrial, lapp_point, ltest, ljacobian);
    else
      localAssembleJacobianVolume(time, ltrial, llin_point, ltest, mat);
  }

  /**
   * @brief      The jacobian volume integral
   *
   * @param[in]  eg    The entity
   * @param[in]  lfsu  The trial local function basis
   * @param[in]  x     The local coefficient vector
   * @param[in]  lfsv  The test local function basis
   * @param      mat   The local matrix
   *
   * @tparam     EG    The entity
   * @tparam     LFSU  The trial local function basis
   * @tparam     X     The local coefficient vector type
   * @tparam     LFSV  The test local function basis
   * @tparam     M     The local matrix
   */
  void localAssembleJacobianVolume(auto time,
                                   const PDELab::Concept::LocalBasis auto& ltrial,
                                   const PDELab::Concept::LocalConstContainer auto& llin_point,
                                   const PDELab::Concept::LocalBasis auto& ltest,
                                   auto& ljacobian) noexcept
  {
    if (ltrial.size() == 0)
      return;

    const auto& entity = ltrial.element();
    const auto& geo = entity.geometry();
    using Geometry = std::decay_t<decltype(geo)>;
    std::optional<typename Geometry::JacobianInverse> geojacinv_opt;

    _local_values_in->time = time;
    _local_values_in->entity_volume = geo.volume();
    _local_values_in->in_volume = 1;

    // Call the gridContext to update the gmsh_id
    _local_values_in->update_gmsh_id(entity);
    // Call the gridContext to update the grid variables
    _local_values_in->update_grid_data(entity);

    auto intorder = integrationOrder(ltrial);
    const auto& quad_rule = QuadratureRules<DF, dim>::rule(geo.type(), intorder);

    // loop over quadrature points
    for (std::size_t q = 0; q != quad_rule.size(); ++q) {
      const auto [position, weight] = quad_rule[q];
      if (not geojacinv_opt or not geo.affine())
        geojacinv_opt.emplace(geo.jacobianInverse(position));
      const auto& geojacinv = *geojacinv_opt;
      _local_values_in->position = geo.global(position);
      auto factor = weight * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      forEachLeafNode(ltrial.tree(), [&](const auto& node) {
        if (node.size() == 0)
          return;
        auto& value = _local_values_in->get_value(node);
        auto& gradient = _local_values_in->get_gradient(node);
        value = 0.;
        gradient = 0.;
        _fe_cache->bind(node.finiteElement(), quad_rule);
        const auto& phi = _fe_cache->evaluateFunction(q);
        const auto& jacphi = _fe_cache->evaluateJacobian(q);
        for (std::size_t dof = 0; dof != node.size(); ++dof) {
            value     += llin_point(node, dof) * phi[dof];
            gradient  += llin_point(node, dof) * (jacphi[dof] * geojacinv)[0];
        }
      });

      // contribution for each component
      forEachLeafNode(ltest.tree(), [&](const auto& ltest_node) {
        if (ltest_node.size() == 0)
          return;
        const auto& eq =
          _local_values_in->get_equation(PDELab::containerEntry(ltrial.tree(), ltest_node.path()));

        _fe_cache->bind(ltest_node.finiteElement(), quad_rule);
        const auto& psi = _fe_cache->evaluateFunction(q);
        const auto& jacpsi = _fe_cache->evaluateJacobian(q);

        // accumulate reaction part into jacobian
        if (eq.reaction) {
          for (const auto& jacobian_entry : eq.reaction.compartment_jacobian) {
            auto jac = jacobian_entry();
            const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
            _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule);
            const auto& phi = _fe_cache->evaluateFunction(q);
            for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
              for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
                ljacobian.accumulate(
                  ltest_node, dof_i, wrt_lbasis, dof_j, -jac * phi[dof_i] * psi[dof_j] * factor);
          }
        }

        if (eq.storage) {
          auto stg = eq.storage();
          const auto& wrt_lbasis = eq.to_local_basis_node(ltrial);
          _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule);
          const auto& phi = _fe_cache->evaluateFunction(q);
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
              ljacobian.accumulate(
                ltest_node, dof_i, wrt_lbasis, dof_j, stg * phi[dof_i] * psi[dof_j] * factor);

          for (const auto& jacobian_entry : eq.storage.compartment_jacobian) {
            auto jac = jacobian_entry();
            const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
            _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule);
            const auto& phi = _fe_cache->evaluateFunction(q);
            for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
              for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
                ljacobian.accumulate(ltest_node,
                                     dof_i,
                                     wrt_lbasis,
                                     dof_j,
                                     jac * eq.value * phi[dof_i] * psi[dof_j] * factor);
          }
        }

        if (eq.velocity) {
          auto vel = eq.velocity();
          const auto& wrt_lbasis = eq.to_local_basis_node(ltrial);
          _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule);
          const auto& phi = _fe_cache->evaluateFunction(q);
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i) {
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
              ljacobian.accumulate(ltest_node,
                                   dof_i,
                                   wrt_lbasis,
                                   dof_j,
                                   -dot(vel * phi[dof_i][0], (jacpsi[dof_j] * geojacinv)[0]) * factor);
          }

          // accumulate jacobian for non-linear terms
          for (const auto& jacobian_entry : eq.velocity.compartment_jacobian) {
            auto adv_flux = jacobian_entry() * eq.value[0];
            const auto& jac_wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
            _fe_cache->bind(jac_wrt_lbasis.finiteElement(), quad_rule);
            const auto& phi = _fe_cache->evaluateFunction(q);
            for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
              for (std::size_t dof_j = 0; dof_j != jac_wrt_lbasis.size(); ++dof_j)
                ljacobian.accumulate(ltest_node,
                                     dof_i,
                                     jac_wrt_lbasis,
                                     dof_j,
                                     -phi[dof_i] * dot(adv_flux, (jacpsi[dof_j] * geojacinv)[0]) * factor);
          }
        }

        // accumulate cross-diffusion part into jacobian
        for (const auto& diffusion : eq.cross_diffusion) {
          const auto& wrt_lbasis = diffusion.wrt.to_local_basis_node(ltrial);
          _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule);
          const auto& jacphi = _fe_cache->evaluateJacobian(q);
          // by product rule
          // accumulate linear term
          for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i) {
            auto diffusive_flux = diffusion((jacphi[dof_i] * geojacinv)[0]);
            for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
              ljacobian.accumulate(
                ltest_node, dof_i, wrt_lbasis, dof_j, dot(diffusive_flux, (jacpsi[dof_j] * geojacinv)[0]) * factor);
          }

          // accumulate jacobian for non-linear terms
          for (const auto& jacobian_entry : diffusion.compartment_jacobian) {
            auto diffusive_flux = jacobian_entry(jacobian_entry.wrt.gradient);
            const auto& jac_wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
            _fe_cache->bind(jac_wrt_lbasis.finiteElement(), quad_rule);
            const auto& phi = _fe_cache->evaluateFunction(q);
            for (std::size_t dof_i = 0; dof_i != ltest_node.size(); ++dof_i)
              for (std::size_t dof_j = 0; dof_j != jac_wrt_lbasis.size(); ++dof_j)
                ljacobian.accumulate(ltest_node,
                                     dof_i,
                                     jac_wrt_lbasis,
                                     dof_j,
                                     phi[dof_i] * dot(diffusive_flux, (jacpsi[dof_j] * geojacinv)[0]) * factor);
          }
        }
      });
    }

    _local_values_in->clear();
    _local_values_in->in_volume = 1;
  }

  template<PDELab::Concept::LocalBasis LocalBasisTrial, PDELab::Concept::LocalBasis LocalBasisTest>
  void localAssembleSkeleton(const Dune::Concept::Intersection auto& intersection,
                             auto time,
                             const LocalBasisTrial& ltrial_in,
                             const PDELab::Concept::LocalConstContainer auto& lcoefficients_in,
                             const LocalBasisTest& ltest_in,
                             const LocalBasisTrial& ltrial_out,
                             const PDELab::Concept::LocalConstContainer auto& lcoefficients_out,
                             const LocalBasisTest& ltest_out,
                             PDELab::Concept::LocalMutableContainer auto& lresidual_in,
                             PDELab::Concept::LocalMutableContainer auto& lresidual_out) noexcept
  {
    if (ltrial_in.size() == 0 and ltrial_out.size() == 0)
      return;

    const auto& entity_i = intersection.inside();
    const auto& geo_i = entity_i.geometry();
    const auto& geo_in_i = intersection.geometryInInside();
    // in case of a boundary, outside objects are an alias of the inside ones
    const auto& entity_o = intersection.neighbor() ? intersection.outside() : entity_i;
    const auto& geo_in_o = intersection.neighbor() ? intersection.geometryInOutside() : geo_in_i;
    const auto& geo_o = entity_o.geometry();

    auto domain_set_i = subDomains(entity_i);
    auto domain_set_o = subDomains(entity_o);

    if (intersection.neighbor() and domain_set_i == domain_set_o)
      return; // not an intersection case

    auto geo_f = intersection.geometry();

    using Geometry = std::decay_t<decltype(entity_i.geometry())>;
    std::optional<typename Geometry::JacobianInverse> geojacinv_opt_i, geojacinv_opt_o;

    _local_values_in->time = _local_values_out->time = time;
    _local_values_in->entity_volume = _local_values_out->entity_volume = geo_f.volume();
    _local_values_in->in_boundary = _local_values_out->in_boundary = static_cast<double>(not intersection.neighbor());
    _local_values_in->in_skeleton = _local_values_out->in_skeleton = static_cast<double>(intersection.neighbor());

    _outflow_i.clear();
    _outflow_o.clear();

    // collect ouflow part for the inside compartment
    if (ltrial_in.size() != 0)
      forEachLeafNode(ltest_in.tree(), [&](const auto& ltest_node_in, auto path) {
        // evaluate outflow only if component exists on this side but not on
        // the outside entity
        const auto& eq =
          _local_values_in->get_equation(PDELab::containerEntry(ltrial_in.tree(), path));
        auto compartment_i = eq.path[Indices::_1];
        if (ltest_node_in.size() == 0 or eq.outflow.empty())
          return;

        // accumulate outflow part into residual
        if (intersection.neighbor()) { // interior skeleton case
          if (domain_set_o.contains(_compartment2domain[compartment_i]))
            return;
          for (std::size_t compartment_o = 0; compartment_o != _compartment2domain.size();
               ++compartment_o) {
            auto domain_o = _compartment2domain[compartment_o];
            if (compartment_i != compartment_o and
                (domain_set_o.contains(domain_o) or domain_set_i.contains(domain_o)) and
                eq.outflow[compartment_o])
              _outflow_i.emplace_back(eq.outflow[compartment_o], eq);
          }
        } else if (eq.outflow[compartment_i]) { // boundary case
          _outflow_i.emplace_back(eq.outflow[compartment_i], eq);
        }
      });

    // collect ouflow part for the outside compartment
    if (intersection.neighbor() and ltrial_out.size() != 0)
      forEachLeafNode(ltest_out.tree(), [&](const auto& ltest_node_out, auto path) {
        // evaluate outflow only if component exists on this side but not on
        // the outside entity
        const auto& eq =
          _local_values_out->get_equation(PDELab::containerEntry(ltrial_out.tree(), path));
        auto compartment_o = eq.path[Indices::_1];
        if (ltest_node_out.size() == 0 or
            domain_set_i.contains(_compartment2domain[compartment_o]) or eq.outflow.empty())
          return;

        // accumulate outflow part into residual (interior skeleton case)
        for (std::size_t compartment_i = 0; compartment_i != _compartment2domain.size();
             ++compartment_i) {
          auto domain_i = _compartment2domain[compartment_i];
          if (compartment_i != compartment_o and
              (domain_set_o.contains(domain_i) or domain_set_i.contains(domain_i)) and
              eq.outflow[compartment_i])
            _outflow_o.emplace_back(eq.outflow[compartment_i], eq);
        }
      });

    if (_outflow_i.empty() and _outflow_o.empty())
      return;

    auto intorder = integrationOrder(ltrial_in, ltrial_out);
    const auto& quad_rule = QuadratureRules<DF, dim-1>::rule(geo_f.type(), intorder);

    // loop over quadrature points
    for (std::size_t q = 0; q != quad_rule.size(); ++q) {
      const auto [position_f, weight] = quad_rule[q];
      auto factor = weight * geo_f.integrationElement(position_f);
      auto normal = intersection.unitOuterNormal(position_f);
      _local_values_in->position = _local_values_out->position = geo_f.global(position_f);
      _local_values_out->normal = -(_local_values_in->normal = normal);

      if (not _outflow_i.empty()) {
        auto quad_proj = [&](auto quad_pos){ return geo_in_i.global(quad_pos); };
        const auto position_i = quad_proj(position_f);
        if (not geojacinv_opt_i or not geo_i.affine())
          geojacinv_opt_i.emplace(geo_i.jacobianInverse(position_i));
        const auto& geojacinv_i = *geojacinv_opt_i;

        // evaluate concentrations at quadrature point (inside part)
        forEachLeafNode(ltrial_in.tree(), [&](const auto& node_in, auto path) {
          const auto& node_out = PDELab::containerEntry(ltrial_out.tree(), path);
          if (node_in.size() == 0 and node_out.size() == 0)
            return;
          // take inside values unless they only exists outside
          const auto& node = (node_in.size() != 0) ? node_in : node_out;
          const auto& lcoefficients = (node_in.size() != 0) ? lcoefficients_in : lcoefficients_out;
          auto& value = _local_values_in->get_value(node);
          auto& gradient = _local_values_in->get_gradient(node);
          value = 0.;
          gradient = 0.;
          _fe_cache->bind(node.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& phi = _fe_cache->evaluateFunction(q);
          const auto& jacphi = _fe_cache->evaluateJacobian(q);
          for (std::size_t dof = 0; dof != node.size(); ++dof) {
            value     += lcoefficients(node, dof) * phi[dof];
            gradient  += lcoefficients(node, dof) * (jacphi[dof] * geojacinv_i)[0];
          }
        });

        // contribution for each component
        for (const auto& [outflow_i, source_i] : _outflow_i) {
          auto outflow = std::invoke(outflow_i);
          const auto& ltest_node_in = source_i.to_local_basis_node(ltest_in);
          _fe_cache->bind(ltest_node_in.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& psi_i = _fe_cache->evaluateFunction(q);
          for (std::size_t dof = 0; dof != ltest_node_in.size(); ++dof)
            lresidual_in.accumulate(ltest_node_in, dof, outflow * psi_i[dof] * factor);
        }
      }

      if (not _outflow_o.empty()) {
        auto quad_proj = [&](auto quad_pos){ return geo_in_o.global(quad_pos); };
        const auto position_o = quad_proj(position_f);
        if (not geojacinv_opt_o or not geo_o.affine())
          geojacinv_opt_o.emplace(geo_o.jacobianInverse(position_o));
        const auto& geojacinv_o = *geojacinv_opt_o;

        // evaluate concentrations at quadrature point (outside part)
        forEachLeafNode(ltrial_out.tree(), [&](const auto& node_out, auto path) {
          const auto& node_in = PDELab::containerEntry(ltrial_in.tree(), path);
          // take outside values unless they only exists inside
          const auto& node = (node_out.size() != 0) ? node_out : node_in;
          const auto& lcoefficients = (node_out.size() != 0) ? lcoefficients_out : lcoefficients_in;
          auto& value = _local_values_out->get_value(node);
          auto& gradient = _local_values_out->get_gradient(node);
          value = 0.;
          gradient = 0.;
          _fe_cache->bind(node.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& phi = _fe_cache->evaluateFunction(q);
          const auto& jacphi = _fe_cache->evaluateJacobian(q);
          for (std::size_t dof = 0; dof != node.size(); ++dof) {
            value     += lcoefficients(node, dof) * phi[dof];
            gradient  += lcoefficients(node, dof) * (jacphi[dof] * geojacinv_o)[0];
          }
        });

        // contribution for each component
        for (const auto& [outflow_o, source_o] : _outflow_o) {
          auto outflow = std::invoke(outflow_o);
          const auto& ltest_node_out = source_o.to_local_basis_node(ltest_out);
          _fe_cache->bind(ltest_node_out.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& psi_o = _fe_cache->evaluateFunction(q);
          for (std::size_t dof = 0; dof != ltest_node_out.size(); ++dof)
            lresidual_out.accumulate(ltest_node_out, dof, outflow * psi_o[dof] * factor);
        }
      }
    }

    _local_values_in->clear();
    _local_values_out->clear();
    _local_values_in->in_boundary = _local_values_in->in_skeleton = 0;
    _local_values_out->in_boundary = _local_values_out->in_skeleton = 0;
  }

  void localAssembleJacobianSkeleton(
    const Dune::Concept::Intersection auto& intersection,
    auto time,
    const PDELab::Concept::LocalBasis auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer auto& llin_point_in,
    const PDELab::Concept::LocalBasis auto& ltest_in,
    const PDELab::Concept::LocalBasis auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer auto& llin_point_out,
    const PDELab::Concept::LocalBasis auto& ltest_out,
    auto& ljacobian_in_in,
    auto& ljacobian_in_out,
    auto& ljacobian_out_in,
    auto& ljacobian_out_out) noexcept
  {
    if (ltrial_in.size() == 0 and ltrial_out.size() == 0)
      return;

    const auto& entity_i = intersection.inside();
    const auto& geo_i = entity_i.geometry();
    const auto& geo_in_i = intersection.geometryInInside();
    // in case of a boundary, outside objects are an alias of the inside ones
    const auto& entity_o = intersection.neighbor() ? intersection.outside() : entity_i;
    const auto& geo_in_o = intersection.neighbor() ? intersection.geometryInOutside() : geo_in_i;
    const auto& geo_o = entity_o.geometry();

    auto domain_set_i = subDomains(entity_i);
    auto domain_set_o = subDomains(entity_o);

    if (intersection.neighbor() and domain_set_i == domain_set_o)
      return; // not an intersection case

    auto geo_f = intersection.geometry();

    using Geometry = std::decay_t<decltype(entity_i.geometry())>;
    std::optional<typename Geometry::JacobianInverse> geojacinv_opt_i, geojacinv_opt_o;

    _local_values_in->time = _local_values_out->time = time;
    _local_values_in->entity_volume = _local_values_out->entity_volume = geo_f.volume();
    _local_values_in->in_boundary = _local_values_out->in_boundary = static_cast<double>(not intersection.neighbor());
    _local_values_in->in_skeleton = _local_values_out->in_skeleton = static_cast<double>(intersection.neighbor());

    _outflow_i.clear();
    _outflow_o.clear();

    // collect ouflow part for the inside compartment
    if (ltrial_in.size() != 0)
      forEachLeafNode(ltest_in.tree(), [&](const auto& ltest_node_in, auto path) {
        // evaluate outflow only if component exists on this side but not on
        // the outside entity
        const auto& eq =
          _local_values_in->get_equation(PDELab::containerEntry(ltrial_in.tree(), path));
        auto compartment_i = eq.path[Indices::_1];
        if (ltest_node_in.size() == 0 or eq.outflow.empty())
          return;

        // accumulate outflow part into residual
        if (intersection.neighbor()) { // interior skeleton case
          if (domain_set_o.contains(_compartment2domain[compartment_i]))
            return;
          for (std::size_t compartment_o = 0; compartment_o != _compartment2domain.size();
               ++compartment_o) {
            auto domain_o = _compartment2domain[compartment_o];
            if (compartment_i != compartment_o and
                (domain_set_o.contains(domain_o) or domain_set_i.contains(domain_o)) and
                eq.outflow[compartment_o])
              if (not eq.outflow[compartment_o].compartment_jacobian.empty())
                _outflow_i.emplace_back(eq.outflow[compartment_o], eq);
          }
        } else if (eq.outflow[compartment_i]) { // boundary case
          if (not eq.outflow[compartment_i].compartment_jacobian.empty())
            _outflow_i.emplace_back(eq.outflow[compartment_i], eq);
        }
      });

    // collect ouflow part for the outside compartment
    if (intersection.neighbor() and ltrial_out.size() != 0)
      forEachLeafNode(ltest_out.tree(), [&](const auto& ltest_node_out, auto path) {
        // evaluate outflow only if component exists on this side but not on
        // the outside entity
        const auto& eq =
          _local_values_out->get_equation(PDELab::containerEntry(ltrial_out.tree(), path));
        auto compartment_o = eq.path[Indices::_1];
        if (ltest_node_out.size() == 0 or
            domain_set_i.contains(_compartment2domain[compartment_o]) or eq.outflow.empty())
          return;

        // accumulate outflow part into residual (interior skeleton case)
        for (std::size_t compartment_i = 0; compartment_i != _compartment2domain.size();
             ++compartment_i) {
          auto domain_i = _compartment2domain[compartment_i];
          if (compartment_i != compartment_o and
              (domain_set_o.contains(domain_i) or domain_set_i.contains(domain_i)) and
              eq.outflow[compartment_i])
            if (not eq.outflow[compartment_i].compartment_jacobian.empty())
              _outflow_o.emplace_back(eq.outflow[compartment_i], eq);
        }
      });

    if (_outflow_i.empty() and _outflow_o.empty())
      return;

    auto intorder = integrationOrder(ltrial_in, ltrial_out);
    const auto& quad_rule = QuadratureRules<DF, dim-1>::rule(geo_f.type(), intorder);

    // loop over quadrature points
    for (std::size_t q = 0; q != quad_rule.size(); ++q) {
      const auto [position_f, weight] = quad_rule[q];
      auto factor = weight * geo_f.integrationElement(position_f);
      auto normal = intersection.unitOuterNormal(position_f);
      _local_values_in->position = _local_values_out->position = geo_f.global(position_f);
      _local_values_out->normal = -(_local_values_in->normal = normal);

      if (not _outflow_i.empty()) {
        auto quad_proj = [&](auto quad_pos){ return geo_in_i.global(quad_pos); };
        const auto position_i = quad_proj(position_f);
        if (not geojacinv_opt_i or not geo_i.affine())
          geojacinv_opt_i.emplace(geo_i.jacobianInverse(position_i));
        const auto& geojacinv_i = *geojacinv_opt_i;

        // evaluate concentrations at quadrature point (inside part)
        forEachLeafNode(ltrial_in.tree(), [&](const auto& node_in, auto path) {
          const auto& node_out = PDELab::containerEntry(ltrial_out.tree(), path);
          if (node_in.size() == 0 and node_out.size() == 0)
            return;
          // take inside values unless they only exists outside
          const auto& node = (node_in.size() != 0) ? node_in : node_out;
          const auto& llin_point = (node_in.size() != 0) ? llin_point_in : llin_point_out;
          auto& value = _local_values_in->get_value(node);
          auto& gradient = _local_values_in->get_gradient(node);
          value = 0.;
          gradient = 0.;
          _fe_cache->bind(node.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& phi = _fe_cache->evaluateFunction(q);
          const auto& jacphi = _fe_cache->evaluateJacobian(q);
          for (std::size_t dof = 0; dof != node.size(); ++dof) {
            value     += llin_point(node, dof) * phi[dof];
            gradient  += llin_point(node, dof) * (jacphi[dof] * geojacinv_i)[0];
          }
        });

        // contribution for each component
        for (const auto& [outflow_i, source_i] : _outflow_i) {
          const auto& ltest_node_in = source_i.to_local_basis_node(ltest_in);
          _fe_cache->bind(ltest_node_in.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& psi = _fe_cache->evaluateFunction(q);
          for (const auto& jacobian_entry : outflow_i.compartment_jacobian) {
            auto jac = jacobian_entry();
            bool do_self_basis = jacobian_entry.wrt.to_local_basis_node(ltrial_out).size() == 0;
            const auto& ltrial = do_self_basis ? ltrial_in : ltrial_out;
            auto& ljacobian = do_self_basis ? ljacobian_in_in : ljacobian_in_out;
            const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
            _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
            const auto& phi = _fe_cache->evaluateFunction(q);
            for (std::size_t dof_i = 0; dof_i != ltest_node_in.size(); ++dof_i)
              for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
                ljacobian.accumulate(ltest_node_in,
                                     dof_i,
                                     wrt_lbasis,
                                     dof_j,
                                     jac * phi[dof_i] * psi[dof_j] * factor);
          }
        }
      }

      if (not _outflow_o.empty()) {
        auto quad_proj = [&](auto quad_pos){ return geo_in_o.global(quad_pos); };
        const auto position_o = quad_proj(position_f);
        if (not geojacinv_opt_o or not geo_o.affine())
          geojacinv_opt_o.emplace(geo_o.jacobianInverse(position_o));
        const auto& geojacinv_o = *geojacinv_opt_o;

        // evaluate concentrations at quadrature point (outside part)
        forEachLeafNode(ltrial_out.tree(), [&](const auto& node_out, auto path) {
          const auto& node_in = PDELab::containerEntry(ltrial_in.tree(), path);
          // take outside values unless they only exists inside
          const auto& node = (node_out.size() != 0) ? node_out : node_in;
          const auto& llin_point = (node_out.size() != 0) ? llin_point_out : llin_point_in;
          auto& value = _local_values_out->get_value(node);
          auto& gradient = _local_values_out->get_gradient(node);
          value = 0.;
          gradient = 0.;
          _fe_cache->bind(node.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& phi = _fe_cache->evaluateFunction(q);
          const auto& jacphi = _fe_cache->evaluateJacobian(q);
          for (std::size_t dof = 0; dof != node.size(); ++dof) {
            value     += llin_point(node, dof) * phi[dof];
            gradient  += llin_point(node, dof) * (jacphi[dof] * geojacinv_o)[0];
          }
        });

        for (const auto& [outflow_o, source_o] : _outflow_o) {
          const auto& ltest_node_out = source_o.to_local_basis_node(ltest_out);
          _fe_cache->bind(ltest_node_out.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
          const auto& psi = _fe_cache->evaluateFunction(q);
          for (const auto& jacobian_entry : outflow_o.compartment_jacobian) {
            auto jac = jacobian_entry();
            bool do_self_basis = jacobian_entry.wrt.to_local_basis_node(ltrial_in).size() == 0;
            const auto& ltrial = do_self_basis ? ltrial_out : ltrial_in;
            auto& ljacobian = do_self_basis ? ljacobian_out_out : ljacobian_out_in;
            const auto& wrt_lbasis = jacobian_entry.wrt.to_local_basis_node(ltrial);
            _fe_cache->bind(wrt_lbasis.finiteElement(), quad_rule, quad_proj, not intersection.conforming());
            const auto& phi = _fe_cache->evaluateFunction(q);
            for (std::size_t dof_i = 0; dof_i != ltest_node_out.size(); ++dof_i)
              for (std::size_t dof_j = 0; dof_j != wrt_lbasis.size(); ++dof_j)
                ljacobian.accumulate(ltest_node_out,
                                     dof_i,
                                     wrt_lbasis,
                                     dof_j,
                                     jac * phi[dof_i] * psi[dof_j] * factor);
          }
        }
      }
    }

    _local_values_in->clear();
    _local_values_out->clear();
    _local_values_in->in_boundary = _local_values_in->in_skeleton = 0;
    _local_values_out->in_boundary = _local_values_out->in_skeleton = 0;
  }

  void localAssembleJacobianSkeletonApply(
    const Dune::Concept::Intersection auto& intersection,
    auto time,
    const PDELab::Concept::LocalBasis auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer auto& lapp_point_in,
    const PDELab::Concept::LocalBasis auto& ltest_in,
    const PDELab::Concept::LocalBasis auto& ltrial_out,
    const PDELab::Concept::LocalConstContainer auto& llin_point_out,
    const PDELab::Concept::LocalConstContainer auto& lapp_point_out,
    const PDELab::Concept::LocalBasis auto& ltest_out,
    PDELab::Concept::LocalMutableContainer auto& ljacobian_in,
    PDELab::Concept::LocalMutableContainer auto& ljacobian_out) noexcept
  {
    PseudoJacobian mat_ii{ ljacobian_in, lapp_point_in };
    PseudoJacobian mat_io{ ljacobian_in, lapp_point_out };
    PseudoJacobian mat_oi{ ljacobian_out, lapp_point_in };
    PseudoJacobian mat_oo{ ljacobian_out, lapp_point_out };
    if (localAssembleIsLinear())
      localAssembleSkeleton(intersection,
                            time,
                            ltrial_in,
                            lapp_point_in,
                            ltest_in,
                            ltrial_out,
                            lapp_point_out,
                            ltest_out,
                            ljacobian_in,
                            ljacobian_out);
    else
      localAssembleJacobianSkeleton(intersection,
                                    time,
                                    ltrial_in,
                                    llin_point_in,
                                    ltest_in,
                                    ltrial_out,
                                    llin_point_out,
                                    ltest_out,
                                    mat_ii,
                                    mat_io,
                                    mat_oi,
                                    mat_oo);
  }

  void localAssembleBoundary(const Dune::Concept::Intersection auto& intersection,
                             auto time,
                             const PDELab::Concept::LocalBasis auto& ltrial_in,
                             const PDELab::Concept::LocalConstContainer auto& lcoefficients_in,
                             const PDELab::Concept::LocalBasis auto& ltest_in,
                             PDELab::Concept::LocalMutableContainer auto& lresidual_in) noexcept
  {
    localAssembleSkeleton(intersection,
                          time,
                          ltrial_in,
                          lcoefficients_in,
                          ltest_in,
                          ltrial_in,
                          lcoefficients_in,
                          ltest_in,
                          lresidual_in,
                          lresidual_in);
  }

  void localAssembleJacobianBoundary(const Dune::Concept::Intersection auto& intersection,
                                     auto time,
                                     const PDELab::Concept::LocalBasis auto& ltrial_in,
                                     const PDELab::Concept::LocalConstContainer auto& llin_point_in,
                                     const PDELab::Concept::LocalBasis auto& ltest_in,
                                     auto& ljacobian_ii) noexcept
  {
    localAssembleJacobianSkeleton(intersection,
                                  time,
                                  ltrial_in,
                                  llin_point_in,
                                  ltest_in,
                                  ltrial_in,
                                  llin_point_in,
                                  ltest_in,
                                  ljacobian_ii,
                                  ljacobian_ii,
                                  ljacobian_ii,
                                  ljacobian_ii);
  }

  void localAssembleJacobianBoundaryApply(
    const Dune::Concept::Intersection auto& intersection,
    auto time,
    const PDELab::Concept::LocalBasis auto& ltrial_in,
    const PDELab::Concept::LocalConstContainer auto& llin_point_in,
    const PDELab::Concept::LocalConstContainer auto& lapp_point_in,
    const PDELab::Concept::LocalBasis auto& ltest_in,
    PDELab::Concept::LocalMutableContainer auto& ljacobian_in) noexcept
  {
    PseudoJacobian mat_ii{ ljacobian_in, lapp_point_in };
    if (localAssembleIsLinear())
      localAssembleBoundary(intersection, time, ltrial_in, lapp_point_in, ltest_in, ljacobian_in);
    else
      localAssembleJacobianBoundary(intersection, time, ltrial_in, llin_point_in, ltest_in, mat_ii);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_CG_HH
