#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_HH

// #include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/quadraturerules.hh>
// #include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
// #include <dune/pdelab/finiteelement/localbasiscache.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/common/fvector.hh>

#include <dune/copasi/pdelab_expression_adapter.hh>

namespace Dune::Copasi {

template<class GridView, class LocalFiniteElement>
class LocalOperatorDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
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
  //! number of total components
  const std::size_t _components;
  //! reference to a quadrature rule
  const QuadratureRule<RF, dim>& _rule;

  //! basis functions at quadrature points
  std::vector<std::vector<Range>> _phihat;
  //! basis function gradients at quadrature points
  std::vector<std::vector<Jacobian>> _gradhat;

  ExpressionToGridFunctionAdapter<GridView, RF> _diffusion_gf;
  mutable ExpressionToGridFunctionAdapter<GridView, RF> _reaction_gf;
  mutable ExpressionToGridFunctionAdapter<GridView, RF> _jacobian_gf;

  Logging::Logger _logger;

  std::set<std::pair<std::size_t, std::size_t>> _component_pattern;

public:
  //! pattern assembly flags
  static constexpr bool doPatternVolume = true;

  //! residual assembly flags
  static constexpr bool doAlphaVolume = true;

  LocalOperatorDiffusionReaction(GridView& grid_view,
                                 const ParameterTree& config,
                                 const LocalFiniteElement& finite_element)
    : _basis_size(finite_element.localBasis().size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _rule(QuadratureRules<RF, dim>::rule(finite_element.type(),
                                           3)) // TODO: make order variable
    , _diffusion_gf(grid_view, config.sub("diffusion"))
    , _reaction_gf(grid_view, config.sub("reaction"), config.sub("reaction"))
    , _jacobian_gf(grid_view,
                   config.sub("reaction.jacobian"),
                   config.sub("reaction"))
    , _logger(Logging::Logging::componentLogger(config, "default"))
  {
    assert(_components == config.sub("diffusion").getValueKeys().size());
    assert(std::pow(_components, 2) ==
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

    auto reaction_config = config.sub("reaction");
    auto jacobian_config = config.sub("reaction.jacobian");

    auto reaction_keys = reaction_config.getValueKeys();
    auto jacobian_keys = jacobian_config.getValueKeys();

    std::sort(reaction_keys.begin(), reaction_keys.end());
    std::sort(jacobian_keys.begin(), jacobian_keys.end());

    std::size_t count = 0;
    for (std::size_t i = 0; i < _components; i++) {
      for (std::size_t j = 0; j < _components; j++, count++) {
        if (i == j) {
          _component_pattern.insert(std::make_pair(i, j));
          continue;
        }
        std::string jacobian =
          jacobian_config.template get<std::string>(jacobian_keys[count]);

        bool do_pattern = true;
        do_pattern &= (jacobian != "0");
        do_pattern &= (jacobian != "0.0");
        do_pattern &= (jacobian != ".0");
        do_pattern &= (jacobian != "0.");
        if (do_pattern)
          _component_pattern.insert(std::make_pair(i, j));
      }
    }

    for (auto i : _component_pattern) {
      _logger.trace("pattern <{},{}>"_fmt, i.first, i.second);
    }

    _logger.debug("LocalOperatorDiffusionReaction constructed"_fmt);
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    auto do_link = [&](std::size_t dof_i, std::size_t dof_j) {
      std::size_t comp_i = dof_i / _basis_size;
      std::size_t comp_j = dof_j / _basis_size;
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };
    for (std::size_t i = 0; i < lfsv.size(); ++i)
      for (std::size_t j = 0; j < lfsu.size(); ++j)
        if (do_link(i, j))
          pattern.addLink(lfsv, i, lfsu, j);
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
    auto x_coeff = [&](const std::size_t& component, const std::size_t& dof) {
      return x(lfsu, component * _basis_size + dof);
    };
    auto z_coeff = [&](const std::size_t& component, const std::size_t& dof) {
      return z(lfsu, component * _basis_size + dof);
    };

    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu, component * _basis_size + dof, value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      const auto& position = _rule[q].position();

      // get diffusion coefficient
      DynamicVector<RF> diffusion;
      _diffusion_gf.evaluate(entity, position, diffusion);

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u(_components);
      for (std::size_t k = 0; k < _components; k++) // loop over components
        for (std::size_t j = 0; j < _basis_size;
             j++) // loop over ansatz functions
          u[k] += x_coeff(k, j) * _phihat[q][j];

      // get diffusion coefficient
      DynamicVector<RF> reaction;
      _reaction_gf.bind(entity, u);
      _reaction_gf.evaluate(entity, position, reaction);

      // compute gradients of basis functions in transformed element
      // (independent of component)
      DynamicVector<FieldVector<RF, dim>> grad(_basis_size);
      for (std::size_t i = 0; i < dim; i++)             // rows of S
        for (std::size_t k = 0; k < dim; k++)           // columns of S
          for (std::size_t j = 0; j < _basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      // contribution for each component
      for (std::size_t k = 0; k < _components; k++) // loop over components
      {
        // compute gradient u_h
        FieldVector<RF, dim> graduh(.0);
        for (std::size_t d = 0; d < dim; d++)           // rows of grad
          for (std::size_t j = 0; j < _basis_size; j++) // columns of grad
            graduh[d] += grad[j][d] * z_coeff(k, j);

        // scalar products
        for (std::size_t d = 0; d < dim; d++) // rows of grad
          for (std::size_t i = 0; i < _basis_size;
               i++) // loop over test functions
            accumulate(k, i, diffusion[k] * grad[i][d] * graduh[d] * factor);

        // reaction term
        for (std::size_t i = 0; i < _basis_size;
             i++) // loop over test functions
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
    auto x_coeff = [&](const std::size_t& component, const std::size_t& dof) {
      return x(lfsu, component * _basis_size + dof);
    };

    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(lfsv,
                     component_i * _basis_size + dof_i,
                     lfsu,
                     component_j * _basis_size + dof_j,
                     value);
    };

    // get entity
    const auto entity = eg.entity();

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      // local stiffness matrix (independent of component)
      std::vector<std::vector<RF>> A(_basis_size);
      std::fill(A.begin(), A.end(), std::vector<RF>(_basis_size));

      const auto& position = _rule[q].position();

      // get diffusion coefficient
      DynamicVector<RF> diffusion;
      _diffusion_gf.evaluate(entity, position, diffusion);

      // get jacobian and determinant
      FieldMatrix<DF, dim, dim> S = geo.jacobianInverseTransposed(position);
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // evaluate concentrations at quadrature point
      DynamicVector<RF> u(_components);
      for (std::size_t k = 0; k < _components; k++) // loop over components
        for (std::size_t j = 0; j < _basis_size;
             j++) // loop over ansatz functions
          u[k] += x_coeff(k, j) * _phihat[q][j];

      // evaluate reaction term
      DynamicVector<RF> jacobian;
      _jacobian_gf.bind(entity, u);
      _jacobian_gf.evaluate(entity, position, jacobian);

      // compute gradients of basis functions in transformed element
      // (independent of component)
      DynamicVector<FieldVector<RF, dim>> grad(_basis_size);
      for (std::size_t i = 0; i < dim; i++)             // rows of S
        for (std::size_t k = 0; k < dim; k++)           // columns of S
          for (std::size_t j = 0; j < _basis_size; j++) // columns of _gradhat
            grad[j][i] += S[i][k] * _gradhat[q][j][0][k];

      // compute grad^T * grad
      for (std::size_t k = 0; k < _components; k++)
        for (std::size_t i = 0; i < _basis_size; i++)
          for (std::size_t d = 0; d < dim; d++)
            for (std::size_t j = 0; j < _basis_size; j++)
              accumulate(
                k, i, k, j, diffusion[k] * grad[i][d] * grad[j][d] * factor);

      int count = 0;
      for (std::size_t k = 0; k < _components; k++)
        for (std::size_t l = 0; l < _components; l++, count++)
          for (std::size_t i = 0; i < _basis_size; i++)
            for (std::size_t j = 0; j < _basis_size; j++)
              accumulate(k,
                         i,
                         l,
                         j,
                         _phihat[q][i] * jacobian[count] * _phihat[q][j] *
                           factor);
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

template<class GridView, class LocalFiniteElement>
class TemporalLocalOperatorDiffusionReaction
  : public Dune::PDELab::LocalOperatorDefaultFlags
  , public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
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
  //! number of total components
  const std::size_t _components;
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
    GridView& grid_view,
    const ParameterTree& config,
    const LocalFiniteElement& finite_element)
    : _basis_size(finite_element.size())
    , _components(config.sub("reaction").getValueKeys().size())
    , _rule(QuadratureRules<RF, dim>::rule(finite_element.type(),
                                           3)) // TODO: make order variable
    , _logger(Logging::Logging::componentLogger(config, "default"))
  {
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

    auto reaction_config = config.sub("reaction");
    auto jacobian_config = config.sub("reaction.jacobian");

    auto reaction_keys = reaction_config.getValueKeys();
    auto jacobian_keys = jacobian_config.getValueKeys();

    std::sort(reaction_keys.begin(), reaction_keys.end());
    std::sort(jacobian_keys.begin(), jacobian_keys.end());

    std::size_t count = 0;
    for (std::size_t i = 0; i < _components; i++) {
      for (std::size_t j = 0; j < _components; j++, count++) {
        if (i == j) {
          _component_pattern.insert(std::make_pair(i, j));
          continue;
        }

        bool do_pattern = true;
        std::string jacobian =
          jacobian_config.template get<std::string>(jacobian_keys[count]);
        do_pattern &= (jacobian.find("0") != std::string::npos);
        do_pattern &= (jacobian.find("0.0") != std::string::npos);
        do_pattern &= (jacobian.find(".0") != std::string::npos);
        do_pattern &= (jacobian.find("0.") != std::string::npos);
        if (do_pattern)
          _component_pattern.insert(std::make_pair(i, j));
      }
    }

    for (auto i : _component_pattern) {
      _logger.trace("pattern <{},{}>"_fmt, i.first, i.second);
    }

    _logger.debug("TemporalLocalOperatorDiffusionReaction constructed"_fmt);
  }

  // define sparsity pattern of operator representation
  template<typename LFSU, typename LFSV, typename LocalPattern>
  void pattern_volume(const LFSU& lfsu,
                      const LFSV& lfsv,
                      LocalPattern& pattern) const
  {
    auto do_link = [&](std::size_t dof_i, std::size_t dof_j) {
      std::size_t comp_i = dof_i / _basis_size;
      std::size_t comp_j = dof_j / _basis_size;
      auto it = _component_pattern.find(std::make_pair(comp_i, comp_j));
      return (it != _component_pattern.end());
    };
    for (std::size_t i = 0; i < lfsv.size(); ++i)
      for (std::size_t j = 0; j < lfsu.size(); ++j)
        if (do_link(i, j))
          pattern.addLink(lfsv, i, lfsu, j);
  }

  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg,
                    const LFSU& lfsu,
                    const X& x,
                    const LFSV& lfsv,
                    R& r) const
  {
    auto x_coeff = [&](const std::size_t& component, const std::size_t& dof) {
      return x(lfsu, component * _basis_size + dof);
    };
    auto accumulate = [&](const std::size_t& component,
                          const std::size_t& dof,
                          const auto& value) {
      r.accumulate(lfsu, component * _basis_size + dof, value);
    };

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      const auto& position = _rule[q].position();
      // get Jacobian and determinant
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // contribution for each component
      for (std::size_t k = 0; k < _components; k++) // loop over components
      {
        // compute value of component
        double u = 0.0;
        for (std::size_t j = 0; j < _basis_size;
             j++) // loop over ansatz functions
          u += x_coeff(k, j) * _phihat[q][j];

        // reaction term
        for (std::size_t i = 0; i < _basis_size;
             i++) // loop over test functions
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
    // assume we receive a power local finite element!
    // auto x_coeff = [&](const std::size_t& component, const std::size_t& dof)
    // {
    //   return x(lfsu, component * _basis_size + dof);
    // };

    auto accumulate = [&](const std::size_t& component_i,
                          const std::size_t& dof_i,
                          const std::size_t& component_j,
                          const std::size_t& dof_j,
                          const auto& value) {
      mat.accumulate(lfsv,
                     component_i * _basis_size + dof_i,
                     lfsu,
                     component_j * _basis_size + dof_j,
                     value);
    };

    // get geometry
    const auto geo = eg.geometry();

    // loop over quadrature points
    for (std::size_t q = 0; q < _rule.size(); q++) {
      const auto& position = _rule[q].position();
      // get Jacobian and determinant
      RF factor = _rule[q].weight() * geo.integrationElement(position);

      // integrate mass matrix
      for (std::size_t k = 0; k < _components; k++)
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