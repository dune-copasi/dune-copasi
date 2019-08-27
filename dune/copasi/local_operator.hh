#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH

#include <dune/copasi/coefficient_mapper.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/pattern.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/common/fvector.hh>

#include <dune/copasi/pdelab_expression_adapter.hh>

namespace Dune::Copasi {

template<class GV, class LFE, class CM = DefaultCoefficientMapper>
class LocalOperatorDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
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

  // coefficient mapper
  mutable CoefficientMapper _coefficient_mapper;
  //! number of basis per component
  const std::size_t _basis_size;
  //! components
  std::vector<std::size_t> _components;
  //! components
  std::vector<std::size_t> _lfs_components;
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
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  template<class T = int,
           class = std::enable_if_t<std::is_same_v<T, int> and
                                      std::is_default_constructible_v<CM>,
                                    int>>
  LocalOperatorDiffusionReaction(GridView grid_view,
                                 const ParameterTree& config,
                                 const LocalFiniteElement& finite_element,
                                 const std::vector<std::size_t>& lsf_components)
    : LocalOperatorDiffusionReaction(grid_view,
                                     config,
                                     finite_element,
                                     lsf_components,
                                     CoefficientMapper{})
  {}

  LocalOperatorDiffusionReaction(GridView grid_view,
                                 const ParameterTree& config,
                                 const LocalFiniteElement& finite_element,
                                 const std::vector<std::size_t>& lsf_components,
                                 const CoefficientMapper& coefficient_mapper)
    : _coefficient_mapper(coefficient_mapper)
    , _basis_size(finite_element.localBasis().size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _lfs_components(lsf_components)
    , _rule(QuadratureRules<RF, dim>::rule(finite_element.type(),
                                           3)) // TODO: make order variable
    , _diffusion_gf(_components.size())
    , _reaction_gf(_components.size())
    , _jacobian_gf(_components.size() * _components.size())
    , _logger(Logging::Logging::componentLogger(config, "model"))
  {
    std::iota(std::next(_components.begin()), _components.end(), 1);

    assert(_components.size() == config.sub("diffusion").getValueKeys().size());
    assert((_components.size() * _components.size()) ==
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

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      const auto& position = _rule[q].position();

      // get diffusion coefficient
      DynamicVector<RF> diffusion(_lfs_components.size());
      for (std::size_t k = 0; k < _lfs_components.size(); k++)
        _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u(_components.size());
      ;
      for (const auto& k : _components)
        for (std::size_t j = 0; j < _basis_size; j++) // ansatz func. loop
          u[k] += _coefficient_mapper(x_coeff_local, k, j) * _phihat[q][j];

      // get diffusion coefficient
      DynamicVector<RF> reaction(_lfs_components.size());
      for (std::size_t k = 0; k < _lfs_components.size(); k++) {
        _reaction_gf[k]->update(u);
        _reaction_gf[k]->evaluate(entity, position, reaction[k]);
      }

      // compute gradients of basis functions in transformed element
      // (independent of component)
      DynamicVector<FieldVector<RF, dim>> grad(_basis_size);
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
        for (std::size_t d = 0; d < dim; d++)           // rows of grad
          for (std::size_t i = 0; i < _basis_size; i++) // test func. loop
            accumulate(k, i, diffusion[k] * grad[i][d] * graduh[d] * factor);

        // reaction term
        for (std::size_t i = 0; i < _basis_size; i++) // test func. loop
          accumulate(k, i, reaction[k] * _phihat[q][i] * factor);
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

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      // local stiffness matrix (independent of component)
      std::vector<std::vector<RF>> A(_basis_size);
      std::fill(A.begin(), A.end(), std::vector<RF>(_basis_size));

      const auto& position = _rule[q].position();

      // get diffusion coefficient
      DynamicVector<RF> diffusion(_lfs_components.size());
      for (std::size_t k = 0; k < _lfs_components.size(); k++)
        _diffusion_gf[k]->evaluate(entity, position, diffusion[k]);

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u(_components.size());
      for (const auto& k : _components)
        for (std::size_t j = 0; j < _basis_size; j++) //  ansatz func. loop
          u[k] += _coefficient_mapper(x_coeff_local, k, j) * _phihat[q][j];

      // evaluate reaction term
      DynamicVector<RF> jacobian(_lfs_components.size() *
                                 _lfs_components.size());
      for (std::size_t k = 0; k < _lfs_components.size(); k++) {
        for (std::size_t l = 0; l < _lfs_components.size(); l++) {
          const auto j = _lfs_components.size() * k + l;
          _jacobian_gf[j]->update(u);
          _jacobian_gf[j]->evaluate(entity, position, jacobian[j]);
        }
      }

      // compute gradients of basis functions in transformed element
      // (independent of component)
      DynamicVector<FieldVector<RF, dim>> grad(_basis_size);
      for (std::size_t i = 0; i < dim; i++)             // rows of S
        for (std::size_t k = 0; k < dim; k++)           // columns of S
          for (std::size_t j = 0; j < _basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      // compute grad^T * grad
      for (std::size_t k = 0; k < _lfs_components.size(); k++)
        for (std::size_t i = 0; i < _basis_size; i++)
          for (std::size_t d = 0; d < dim; d++)
            for (std::size_t j = 0; j < _basis_size; j++)
              accumulate(
                k, i, k, j, diffusion[k] * grad[i][d] * grad[j][d] * factor);

      for (std::size_t k = 0; k < _lfs_components.size(); k++)
        for (std::size_t l = 0; l < _lfs_components.size(); l++) {
          const auto j = _lfs_components.size() * k + l;
          for (std::size_t m = 0; m < _basis_size; m++)
            for (std::size_t n = 0; n < _basis_size; n++)
              accumulate(k,
                         m,
                         l,
                         n,
                         _phihat[q][m] * jacobian[j] * _phihat[q][n] * factor);
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
    jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    jacobian_apply_volume(eg, lfsu, x, x, lfsv, r);
  }
};

template<class GV, class LFE, class CM = DefaultCoefficientMapper>
class TemporalLocalOperatorDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::FullVolumePattern
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
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

  //! world dimension
  static constexpr int dim = LocalBasis::Traits::dimDomain;

  //! range dimension
  static constexpr int dim_range = LocalBasis::Traits::dimRange;

  // this operator only support scalar ranges
  static_assert(dim_range == 1);

  //! coefficient mapper
  mutable CoefficientMapper _coefficient_mapper;
  //! number of basis per component
  const std::size_t _basis_size;
  //! components
  std::vector<std::size_t> _components;
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

  template<class T = int,
           class = std::enable_if_t<std::is_same_v<T, int> and
                                      std::is_default_constructible_v<CM>,
                                    int>>
  TemporalLocalOperatorDiffusionReaction(
    GridView grid_view,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element,
    const std::vector<std::size_t>& lsf_components)
    : TemporalLocalOperatorDiffusionReaction(grid_view,
                                             config,
                                             finite_element,
                                             lsf_components,
                                             CoefficientMapper{})
  {}

  TemporalLocalOperatorDiffusionReaction(
    GridView grid_view,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element,
    const std::vector<std::size_t>& lfs_components,
    const CoefficientMapper& coefficient_mapper)
    : _coefficient_mapper(coefficient_mapper)
    , _basis_size(finite_element.size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _lfs_components(lfs_components)
    , _rule(QuadratureRules<RF, dim>::rule(finite_element.type(),
                                           3)) // TODO: make order variable
    , _logger(Logging::Logging::componentLogger(config, "model"))
  {
    assert(_components.size() == config.sub("diffusion").getValueKeys().size());

    std::iota(std::next(_components.begin()), _components.end(), 1);

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

  template<class States>
  void update(const States& states)
  {
    _coefficient_mapper.update(states);
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

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    _coefficient_mapper.bind(entity);

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
          u += _coefficient_mapper(x_coeff_local, k, j) * _phihat[q][j];

        // reaction term
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
    alpha_volume(eg, lfsu, z, lfsv, r);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void jacobian_apply_volume(const EG& eg,
                             const LFSU& lfsu,
                             const X& x,
                             const LFSV& lfsv,
                             R& r) const
  {
    alpha_volume(eg, lfsu, x, lfsv, r);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH