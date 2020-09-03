#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH

#include <dune/copasi/common/pdelab_expression_adapter.hh>

#include <dune/pdelab/common/quadraturerules.hh>

#include <set>

namespace Dune::Copasi {

/**
 *
 * @tparam     GV    GridView
 * @tparam     LBT    Local basis traits
 */
template<class GV, class LBT>
struct LocalOperatorDiffusionReactionBase
{
  //! grid view
  using GridView = GV;

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

  std::vector<std::shared_ptr<ExpressionAdapter>> _diffusion_gf;
  mutable std::vector<std::shared_ptr<ExpressionAdapter>> _reaction_gf;
  mutable std::vector<std::shared_ptr<ExpressionAdapter>> _jacobian_gf;

  Logging::Logger _logger;

  std::set<std::pair<std::size_t, std::size_t>> _component_pattern;

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
   */
  LocalOperatorDiffusionReactionBase(const ParameterTree& config)
    : _components(config.sub("reaction").getValueKeys().size())
    , _logger(Logging::Logging::componentLogger({}, "model"))
  {
    assert(_components == config.sub("diffusion").getValueKeys().size());
    assert((_components * _components) ==
           config.sub("reaction.jacobian").getValueKeys().size());
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
   */
  LocalOperatorDiffusionReactionBase(const ParameterTree& config,
                                  const GridView& grid_view)
    : LocalOperatorDiffusionReactionBase(config)
  {
    create_pattern_and_gf_expressions(grid_view,config);
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

    for (std::size_t k = 0; k < _components; k++) {
      std::string var = reaction_keys[k];

      std::string d_eq = diffusion_config.template get<std::string>(var);
      std::string r_eq = reaction_config.template get<std::string>(var);

      _diffusion_gf[k] =
        std::make_shared<ExpressionAdapter>(grid_view, d_eq);
      _reaction_gf[k] =
        std::make_shared<ExpressionAdapter>(grid_view, r_eq, true, reaction_keys);

      for (std::size_t l = 0; l < _components; l++) {
        const auto j = _components * k + l;
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
