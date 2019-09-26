#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH

#include <dune/copasi/coefficient_mapper.hh>
#include <dune/copasi/enum.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
#include <dune/pdelab/localoperator/pattern.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/common/fvector.hh>

#include <dune/copasi/pdelab_expression_adapter.hh>

#include <set>

namespace Dune::Copasi {

template<class GV,
         class LFE,
         class CM = DefaultCoefficientMapper,
         JacobianMethod JM = JacobianMethod::Analytical>
class LocalOperatorDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public PDELab::NumericalJacobianVolume<
      LocalOperatorDiffusionReaction<GV, LFE, CM, JM>>
  , public PDELab::NumericalJacobianApplyVolume<
      LocalOperatorDiffusionReaction<GV, LFE, CM, JM>>
{
  //! grid view
  using GridView = GV;

  //! local finite element
  using LocalFiniteElement = LFE;

  //! coefficient mapper
  using CoefficientMapper = CM;

  //! local basis
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  //! domain field
  using DF = typename LocalBasis::Traits::DomainFieldType;

  //! coordinates type
  using Domain = typename LocalBasis::Traits::DomainType;

  //! range field
  using RF = typename LocalBasis::Traits::RangeFieldType;

  //! range type (for the local finite element)
  using Range = typename LocalBasis::Traits::RangeType;

  //! jacobian tpye
  using Jacobian = typename LocalBasis::Traits::JacobianType;

  //! Adapter for dynamic expressions
  using ExpressionAdapter = ExpressionToGridFunctionAdapter<GridView, RF>;

  //! world dimension
  static constexpr int dim = LocalBasis::Traits::dimDomain;

  //! range dimension
  static constexpr int dim_range = LocalBasis::Traits::dimRange;

  // this operator only support scalar ranges
  static_assert(dim_range == 1);

  //! number of basis per component
  const std::size_t _basis_size;
  //! components
  std::size_t _components;
  //! reference to a quadrature rule
  const QuadratureRule<RF, dim>& _rule;

  //! basis functions at quadrature points
  std::vector<std::vector<Range>> _phihat;
  //! basis function gradients at quadrature points
  std::vector<std::vector<Jacobian>> _gradhat;

  std::vector<std::shared_ptr<ExpressionAdapter>> _diffusion_gf;
  mutable std::vector<std::shared_ptr<ExpressionAdapter>> _reaction_gf;
  mutable std::vector<std::shared_ptr<ExpressionAdapter>> _jacobian_gf;

  Logging::Logger _logger;

  std::set<std::pair<std::size_t, std::size_t>> _component_pattern;

public:
  //! components
  std::vector<std::size_t> _lfs_components;

  // coefficient mapper
  mutable CoefficientMapper _coefficient_mapper;

  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  LocalOperatorDiffusionReaction(GridView grid_view,
                                 const ParameterTree& config,
                                 const LocalFiniteElement& finite_element,
                                 std::size_t id_operator)
    : _basis_size(finite_element.localBasis().size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _rule(QuadratureRules<RF, dim>::rule(finite_element.type(),
                                           3)) // TODO: make order variable
    , _diffusion_gf(_components)
    , _reaction_gf(_components)
    , _jacobian_gf(_components * _components)
    , _logger(Logging::Logging::componentLogger(config, "model"))
    , _coefficient_mapper(config.sub("operator"), id_operator)
  {

    auto& config_operator = config.sub("operator");
    auto comp_names = config_operator.getValueKeys();
    std::sort(comp_names.begin(), comp_names.end());

    for (std::size_t j = 0; j < comp_names.size(); j++) {
      std::size_t k = config_operator.template get<std::size_t>(comp_names[j]);
      if (k == id_operator)
        _lfs_components.push_back(j);
    }

    assert(_components == config.sub("diffusion").getValueKeys().size());
    assert((_components * _components) ==
           config.sub("reaction.jacobian").getValueKeys().size());

    _logger.trace("cache finite element evaluations on reference element"_fmt);
    std::vector<Range> phi(_basis_size);
    std::vector<Jacobian> jac(_basis_size);
    for (const auto& point : _rule) {
      const auto& position = point.position();
      _logger.trace("position: {}"_fmt, position);

      const auto& local_basis = finite_element.localBasis();
      local_basis.evaluateFunction(position, phi);
      local_basis.evaluateJacobian(position, jac);

      for (std::size_t i = 0; i < _basis_size; ++i) {
        _logger.trace(" value[{}]: {}"_fmt, i, phi[i]);
        _logger.trace(" jacobian[{}]: {}"_fmt, i, jac[i]);
      }
      _phihat.push_back(phi);
      phi.clear();
      _gradhat.push_back(jac);
      jac.clear();
    }

    auto diffusion_config = config.sub("diffusion");
    auto reaction_config = config.sub("reaction");
    auto jacobian_config = config.sub("reaction.jacobian");

    auto diffusion_keys = diffusion_config.getValueKeys();
    auto reaction_keys = reaction_config.getValueKeys();
    auto jacobian_keys = jacobian_config.getValueKeys();

    std::sort(diffusion_keys.begin(), diffusion_keys.end());
    std::sort(reaction_keys.begin(), reaction_keys.end());
    std::sort(jacobian_keys.begin(), jacobian_keys.end());

    assert(diffusion_keys.size() == reaction_keys.size());
    for (size_t i = 0; i < diffusion_keys.size(); i++)
      assert(diffusion_keys[i] == reaction_keys[i]);

    for (std::size_t k = 0; k < _lfs_components.size(); k++) {
      std::string var = reaction_keys[_lfs_components[k]];

      std::string d_eq = diffusion_config.template get<std::string>(var);
      std::string r_eq = reaction_config.template get<std::string>(var);

      _diffusion_gf[k] =
        std::make_shared<ExpressionToGridFunctionAdapter<GridView, RF>>(
          grid_view, d_eq);
      _reaction_gf[k] =
        std::make_shared<ExpressionToGridFunctionAdapter<GridView, RF>>(
          grid_view, r_eq, reaction_keys);

      for (std::size_t l = 0; l < _lfs_components.size(); l++) {
        const auto j = _lfs_components.size() * k + l;
        std::string j_eq =
          jacobian_config.template get<std::string>(jacobian_keys[j]);

        _jacobian_gf[j] =
          std::make_shared<ExpressionToGridFunctionAdapter<GridView, RF>>(
            grid_view, j_eq, reaction_keys);
        if (k == l) {
          _component_pattern.insert(std::make_pair(k, l));
          continue;
        }

        bool do_pattern = true;
        do_pattern &= (j_eq != "0");
        do_pattern &= (j_eq != "0.0");
        do_pattern &= (j_eq != ".0");
        do_pattern &= (j_eq != "0.");
        if (do_pattern)
          _component_pattern.insert(std::make_pair(k, l));
      }
    }

    for (auto i : _component_pattern) {
      _logger.trace("pattern <{},{}>"_fmt, i.first, i.second);
    }

    _logger.debug("LocalOperatorDiffusionReaction constructed"_fmt);
  }

  template<class States>
  void update(const States& states)
  {
    _coefficient_mapper.update(states);
  }

  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };
    for (std::size_t i = 0; i < lfsv.degree(); ++i)
      for (std::size_t j = 0; j < lfsu.degree(); ++j)
        if (do_link(i, j))
          for (std::size_t k = 0; k < lfsv.child(i).size(); ++k)
            for (std::size_t l = 0; l < lfsu.child(j).size(); ++l)
              pattern.addLink(lfsv.child(i), k, lfsu.child(j), l);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<LocalOperatorDiffusionReaction>::
        jacobian_apply_volume(eg, lfsu, x, z, lfsv, r);
    } else {
      _jacobian_apply_volume(eg, lfsu, x, z, lfsv, r);
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void _jacobian_apply_volume(const EG& eg,
                              const LFSU& lfsu,
                              const X& x,
                              const X& z,
                              const LFSV& lfsv,
                              R& r) const
  {
    // assume we receive a power local finite element!
    auto x_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return x(lfsu.child(component), dof);
    };
    auto z_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return z(lfsu.child(component), dof);
    };

    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu.child(component), dof, value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    _coefficient_mapper.bind(entity);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> diffusion(_lfs_components.size());
    DynamicVector<RF> reaction(_lfs_components.size());
    DynamicVector<FieldVector<RF, dim>> grad(_basis_size);

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {

      const auto& position = _rule[q].position();

      std::fill(u.begin(), u.end(), 0.);
      std::fill(diffusion.begin(), diffusion.end(), 0.);
      std::fill(reaction.begin(), reaction.end(), 0.);
      std::fill(grad.begin(), grad.end(), FieldVector<RF, dim>(0.));

      // get diffusion coefficient
      for (std::size_t k = 0; k < _lfs_components.size(); k++)
        _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      for (std::size_t k = 0; k < _components; k++)
        for (std::size_t j = 0; j < _basis_size; j++) // ansatz func. loop
          u[k] += _coefficient_mapper(x_coeff_local, k, j) * _phihat[q][j];

      // get reaction term
      for (std::size_t k = 0; k < _lfs_components.size(); k++) {
        _reaction_gf[k]->update(u);
        _reaction_gf[k]->evaluate(entity, position, reaction[k]);
      }

      // compute gradients of basis functions in transformed element
      // (independent of component)
      for (std::size_t i = 0; i < dim; i++)             // rows of S
        for (std::size_t k = 0; k < dim; k++)           // columns of S
          for (std::size_t j = 0; j < _basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      // contribution for each component
      for (std::size_t k = 0; k < _lfs_components.size(); k++) {
        // compute gradient u_h
        FieldVector<RF, dim> graduh(.0);
        for (std::size_t d = 0; d < dim; d++)           // rows of grad
          for (std::size_t j = 0; j < _basis_size; j++) // columns of grad
            graduh[d] += grad[j][d] * z_coeff_local(k, j);

        // scalar products
        for (std::size_t i = 0; i < _basis_size; i++) // test func. loop
        {
          typename R::value_type rhs = -reaction[k] * _phihat[q][i];
          for (std::size_t d = 0; d < dim; d++) // rows of grad
            rhs += diffusion[k] * grad[i][d] * graduh[d];
          accumulate(k, i, rhs * factor);
        }
      }
    }
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       M& mat) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianVolume<LocalOperatorDiffusionReaction>::
        jacobian_volume(eg, lfsu, x, lfsv, mat);
      return;
    }
    // assume we receive a power local finite element!
    auto x_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return x(lfsu.child(component), dof);
    };

    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(
        lfsv.child(component_i), dof_i, lfsu.child(component_j), dof_j, value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    _coefficient_mapper.bind(entity);

    DynamicVector<RF> u(_components);
    DynamicVector<RF> diffusion(_lfs_components.size());
    DynamicVector<RF> jacobian(_lfs_components.size() * _lfs_components.size());
    DynamicVector<FieldVector<RF, dim>> grad(_basis_size);

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {

      const auto& position = _rule[q].position();

      std::fill(u.begin(), u.end(), 0.);
      std::fill(diffusion.begin(), diffusion.end(), 0.);
      std::fill(jacobian.begin(), jacobian.end(), 0.);
      std::fill(grad.begin(), grad.end(), FieldVector<RF, dim>(0.));

      // get diffusion coefficient
      for (std::size_t k = 0; k < _lfs_components.size(); k++)
        _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

      // evaluate concentrations at quadrature point
      for (std::size_t k = 0; k < _components; k++)
        for (std::size_t j = 0; j < _basis_size; j++) //  ansatz func. loop
          u[k] += _coefficient_mapper(x_coeff_local, k, j) * _phihat[q][j];

      // evaluate reaction term
      for (std::size_t k = 0; k < _lfs_components.size(); k++) {
        for (std::size_t l = 0; l < _lfs_components.size(); l++) {
          const auto j = _lfs_components.size() * k + l;
          _jacobian_gf[j]->update(u);
          _jacobian_gf[j]->evaluate(entity, position, jacobian[j]);
        }
      }

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // compute gradients of basis functions in transformed element
      // (independent of component)
      for (std::size_t i = 0; i < dim; i++)             // rows of S
        for (std::size_t k = 0; k < dim; k++)           // columns of S
          for (std::size_t j = 0; j < _basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      auto do_link = [&](std::size_t comp_i, std::size_t comp_j) {
        auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
        return (it != _component_pattern.end());
      };

      // compute grad^T * grad
      for (std::size_t k = 0; k < _lfs_components.size(); k++) {
        for (std::size_t l = 0; l < _lfs_components.size(); l++) {
          if (not do_link(k, l))
            continue;
          const auto j = _lfs_components.size() * k + l;
          for (std::size_t m = 0; m < _basis_size; m++) {
            for (std::size_t n = 0; n < _basis_size; n++) {
              typename M::value_type jac =
                -jacobian[j] * _phihat[q][m] * _phihat[q][n];
              if (l == k)
                for (std::size_t d = 0; d < dim; d++)
                  jac += diffusion[k] * grad[m][d] * grad[n][d];
              accumulate(k, m, l, n, jac * factor);
            }
          }
        }
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    _jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<LocalOperatorDiffusionReaction>::
        jacobian_apply_volume(eg, lfsu, x, lfsv, r);
      return;
    }
    _jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }
};

template<class GV, class LFE, JacobianMethod JM = JacobianMethod::Analytical>
class TemporalLocalOperatorDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
  , public Dune::PDELab::NumericalJacobianVolume<
      TemporalLocalOperatorDiffusionReaction<GV, LFE, JM>>
  , public Dune::PDELab::NumericalJacobianApplyVolume<
      TemporalLocalOperatorDiffusionReaction<GV, LFE, JM>>
{
  //! grid view
  using GridView = GV;

  //! local finite element
  using LocalFiniteElement = LFE;

  //! local basis
  using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

  //! domain field
  using DF = typename LocalBasis::Traits::DomainFieldType;

  //! coordinates type
  using Domain = typename LocalBasis::Traits::DomainType;

  //! range field
  using RF = typename LocalBasis::Traits::RangeFieldType;

  //! range type (for the local finite element)
  using Range = typename LocalBasis::Traits::RangeType;

  //! jacobian tpye
  using Jacobian = typename LocalBasis::Traits::JacobianType;

  //! world dimension
  static constexpr int dim = LocalBasis::Traits::dimDomain;

  //! range dimension
  static constexpr int dim_range = LocalBasis::Traits::dimRange;

  // this operator only support scalar ranges
  static_assert(dim_range == 1);

  //! number of basis per component
  const std::size_t _basis_size;
  //! components
  std::size_t _components;
  //! components
  std::vector<std::size_t> _lfs_components;
  //! reference to a quadrature rule
  const QuadratureRule<RF, dim>& _rule;

  //! basis functions at quadrature points
  std::vector<std::vector<Range>> _phihat;

  Logging::Logger _logger;

  std::set<std::pair<std::size_t, std::size_t>> _component_pattern;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  TemporalLocalOperatorDiffusionReaction(
    GridView grid_view,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element,
    std::size_t id_operator)
    : _basis_size(finite_element.size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _rule(QuadratureRules<RF, dim>::rule(finite_element.type(),
                                           3)) // TODO: make order variable
    , _logger(Logging::Logging::componentLogger(config, "model"))
  {
    auto& config_operator = config.sub("operator");
    auto comp_names = config_operator.getValueKeys();
    std::sort(comp_names.begin(), comp_names.end());

    for (std::size_t j = 0; j < comp_names.size(); j++) {
      std::size_t k = config_operator.template get<std::size_t>(comp_names[j]);
      if (k == id_operator)
        _lfs_components.push_back(j);
    }

    assert(_components == config.sub("diffusion").getValueKeys().size());

    _logger.trace("cache finite element evaluations on reference element"_fmt);
    std::vector<Range> phi(_basis_size);
    for (const auto& point : _rule) {
      const auto& position = point.position();
      _logger.trace("position: {}"_fmt, position);

      const auto& local_basis = finite_element.localBasis();
      local_basis.evaluateFunction(position, phi);

      for (std::size_t i = 0; i < _basis_size; ++i)
        _logger.trace(" value[{}]: {}"_fmt, i, phi[i]);

      _phihat.push_back(phi);
      phi.clear();
    }

    _logger.debug("TemporalLocalOperatorDiffusionReaction constructed"_fmt);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    auto x_coeff_local = [&](const std::size_t& component,
                             const std::size_t& dof) {
      return x(lfsu, component * _basis_size + dof);
    };
    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu.child(component), dof, value);
    };

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      const auto& position = _rule[q].position();
      // get Jacobian and determinant
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // contribution for each component
      for (std::size_t k = 0; k < _lfs_components.size();
           k++) // loop over components
      {
        // compute value of component
        double u = 0.0;
        for (std::size_t j = 0; j < _basis_size; j++) // ansatz func. loop
          u += x_coeff_local(k, j) * _phihat[q][j];

        for (std::size_t i = 0; i < _basis_size; i++) // test func. loop
          accumulate(k, i, u * _phihat[q][i] * factor);
      }
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename Mat>
  void jacobian_volume(const EG& eg,
                       const LFSU& lfsu,
                       const X& x,
                       const LFSV& lfsv,
                       Mat& mat) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianVolume<TemporalLocalOperatorDiffusionReaction>::
        jacobian_volume(eg, lfsu, x, lfsv, mat);
      return;
    }

    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(
        lfsv.child(component_i), dof_i, lfsu.child(component_j), dof_j, value);
    };

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      const auto& position = _rule[q].position();
      // get Jacobian and determinant
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // integrate mass matrix
      for (std::size_t k = 0; k < _lfs_components.size();
           k++) // loop over components
        for (std::size_t i = 0; i < _basis_size; i++)
          for (std::size_t j = 0; j < _basis_size; j++)
            accumulate(k, i, k, j, _phihat[q][i] * _phihat[q][j] * factor);
    }
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const X& z,
                             const LFSV& lfsv,
                             R& r) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<
        TemporalLocalOperatorDiffusionReaction>::jacobian_apply_volume(eg,
                                                                       lfsu,
                                                                       x,
                                                                       z,
                                                                       lfsv,
                                                                       r);
      return;
    }
    alpha_volume(eg, lfsu, z, lfsv, r);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    if constexpr (JM == JacobianMethod::Numerical) {
      PDELab::NumericalJacobianApplyVolume<
        TemporalLocalOperatorDiffusionReaction>::jacobian_apply_volume(eg,
                                                                       lfsu,
                                                                       x,
                                                                       lfsv,
                                                                       r);
      return;
    }
    alpha_volume(eg, lfsu, x, lfsv, r);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH