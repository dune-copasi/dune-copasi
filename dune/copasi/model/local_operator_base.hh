#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH

#include <dune/copasi/common/coefficient_mapper.hh>
#include <dune/copasi/common/pdelab_expression_adapter.hh>

#include <dune/pdelab/common/quadraturerules.hh>

#include <set>

namespace Dune::Copasi {

/**
 *
 * @tparam     GV    GridView
 * @tparam     CM    Local basis traits
 * @tparam     CM    Coefficient mapper
 */
template<class GV, class LBT, class CM = DefaultCoefficientMapper>
struct LocalOperatorDiffusionReactionBase
{
  //! grid view
  using GridView = GV;

  //! coefficient mapper
  using CoefficientMapper = CM;

  //! local basis
  using LocalBasisTraits = LBT;

  //! domain field
  using DF = typename LocalBasisTraits::DomainFieldType;

  //! coordinates type
  using Domain = typename LocalBasisTraits::DomainType;

  //! range field
  using RF = typename LocalBasisTraits::RangeFieldType;

  //! range type (for the local finite element)
  using Range = typename LocalBasisTraits::RangeType;

  //! jacobian tpye
  using Jacobian = typename LocalBasisTraits::JacobianType;

  //! Adapter for dynamic expressions
  using ExpressionAdapter = ExpressionToGridFunctionAdapter<GridView, RF>;

  //! world dimension
  static constexpr int dim = LocalBasisTraits::dimDomain;

  //! range dimension
  static constexpr int dim_range = LocalBasisTraits::dimRange;

  // this operator only support scalar ranges
  static_assert(dim_range == 1);

  //! components
  std::size_t _components;

  bool fe_cached = false;

  //! basis functions at quadrature points
  std::vector<std::vector<Range>> _phihat;
  //! basis function gradients at quadrature points
  std::vector<std::vector<Jacobian>> _gradhat;

  std::vector<std::shared_ptr<ExpressionAdapter>> _diffusion_gf;
  mutable std::vector<std::shared_ptr<ExpressionAdapter>> _reaction_gf;
  mutable std::vector<std::shared_ptr<ExpressionAdapter>> _jacobian_gf;

  Logging::Logger _logger;

  std::set<std::pair<std::size_t, std::size_t>> _component_pattern;

  //! components
  std::vector<std::size_t> _lfs_components;

  // coefficient mapper
  mutable CoefficientMapper _coefficient_mapper_i;
  mutable CoefficientMapper _coefficient_mapper_o;

  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  grid_view       The grid view where this local operator is
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   * @param[in]  rule_order      order of the quadrature rule
   * @param[in]  id_operator     The index of this operator
   */
  LocalOperatorDiffusionReactionBase(const ParameterTree& config,
                                   std::size_t id_operator)
    : _components(config.sub("reaction").getValueKeys().size())
    , _logger(Logging::Logging::componentLogger(config, "model"))
    , _coefficient_mapper_i(config.sub("operator"), id_operator)
    , _coefficient_mapper_o(config.sub("operator"), id_operator)
  {
    assert(_components == config.sub("diffusion").getValueKeys().size());
    assert((_components * _components) ==
           config.sub("reaction.jacobian").getValueKeys().size());

    setup_lfs_components(config,id_operator);
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  grid_view       The grid view where this local operator is
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   * @param[in]  rule_order      order of the quadrature rule
   * @param[in]  id_operator     The index of this operator
   */
  LocalOperatorDiffusionReactionBase(const ParameterTree& config,
                                   std::size_t id_operator,
                                  const GridView& grid_view)
    : LocalOperatorDiffusionReactionBase(config,id_operator)
  {
    create_pattern_and_gf_expressions(grid_view,config);
  }

  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  grid_view       The grid view where this local operator is
   *                             valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   * @param[in]  rule_order      order of the quadrature rule
   * @param[in]  id_operator     The index of this operator
   */
  template<class LocalFiniteElement>
  LocalOperatorDiffusionReactionBase(const ParameterTree& config,
                                   std::size_t id_operator,
                                   const GridView& grid_view,
                                   const LocalFiniteElement& finite_element,
                                   std::size_t rule_order)
    : LocalOperatorDiffusionReactionBase(config,id_operator,grid_view)
  {

    cache_local_finite_element(finite_element,rule_order);
  }

  void setup_lfs_components(const ParameterTree& config, std::size_t id_operator)
  {
    _logger.trace("setup lfs components"_fmt);
    const auto& config_operator = config.sub("operator");
    auto comp_names = config_operator.getValueKeys();
    std::sort(comp_names.begin(), comp_names.end());

    for (std::size_t j = 0; j < comp_names.size(); j++) {
      std::size_t k = config_operator.template get<std::size_t>(comp_names[j]);
      if (k == id_operator)
        _lfs_components.push_back(j);
    }
  }

  template<class LocalFiniteElement>
  void cache_local_finite_element(const LocalFiniteElement& finite_element, std::size_t rule_order)
  {
    _logger.trace("caching local finite element evaluations on reference element"_fmt);

    //! local basis
    using LocalBasis = typename LocalFiniteElement::Traits::LocalBasisType;

    static_assert(std::is_same_v<typename LocalBasis::Traits,LocalBasisTraits>);

    const auto& rule = QuadratureRules<RF, dim>::rule(finite_element.type(),rule_order);

    _phihat.clear();
    _gradhat.clear();
    _phihat.reserve(rule.size());
    _gradhat.reserve(rule.size());

    std::vector<Range> phi(finite_element.localBasis().size());
    std::vector<Jacobian> jac(finite_element.localBasis().size());

    for (const auto& point : rule) {
      const auto& position = point.position();
      _logger.trace("position: {}"_fmt, position);

      const auto& local_basis = finite_element.localBasis();
      local_basis.evaluateFunction(position, phi);
      local_basis.evaluateJacobian(position, jac);

      for (std::size_t i = 0; i < finite_element.localBasis().size(); ++i) {
        _logger.trace(" value[{}]: {}"_fmt, i, phi[i]);
        _logger.trace(" jacobian[{}]: {}"_fmt, i, jac[i]);
      }
      _phihat.push_back(phi);
      phi.clear();
      _gradhat.push_back(jac);
      jac.clear();
    }
  }

  void create_pattern_and_gf_expressions(const GridView& grid_view, const ParameterTree& config)
  {
    _logger.trace("creating pattern and grid function expressions"_fmt);

    _diffusion_gf.resize(_components);
    _reaction_gf.resize(_components);
    _jacobian_gf.resize(_components * _components);

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
        std::make_shared<ExpressionAdapter>(grid_view, d_eq);
      _reaction_gf[k] =
        std::make_shared<ExpressionAdapter>(grid_view, r_eq, true, reaction_keys);

      for (std::size_t l = 0; l < _lfs_components.size(); l++) {
        const auto j = _lfs_components.size() * k + l;
        std::string j_eq =
          jacobian_config.template get<std::string>(jacobian_keys[j]);

        _jacobian_gf[j] =
          std::make_shared<ExpressionAdapter>(grid_view, j_eq, true, reaction_keys);
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
  }

  /**
   * @brief      Updates the coefficient mapper with given states.
   *
   * @param[in]  states  A map from operator index to states
   *
   * @tparam     States  Map from index to states
   */
  template<class States>
  void update(const States& states)
  {
    _coefficient_mapper_i.update(states);
    _coefficient_mapper_o.update(states);
  }

  /**
   * @brief      Sets the time.
   *
   * @param[in]  t     The new time
   */
  void setTime(double t)
  {
    for (auto& gf : _jacobian_gf)
      if (gf)
        gf->set_time(t);
    for (auto& gf : _reaction_gf)
      if (gf)
        gf->set_time(t);
    for (auto& gf : _diffusion_gf)
      if (gf)
        gf->set_time(t);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH