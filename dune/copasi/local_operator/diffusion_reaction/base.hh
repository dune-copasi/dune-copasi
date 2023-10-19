#ifndef DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH
#define DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH


#include <dune/copasi/common/parser_to_grid_function.hh>

#include <dune/assembler/common/trace.hh>
#include <dune/assembler/common/multiindex.hh>
#include <dune/assembler/concepts/discrete_function_space.hh>

#include <dune/pdelab/common/quadraturerules.hh>

#include <set>

namespace Dune::Copasi {

/**
 *
 * @tparam     ES     EntitySet
 * @tparam     LBT    Local basis traits
 */
template<Assembler::Concept::DiscreteFunctionSpace Space, class LBT>
struct LocalOperatorDiffusionReactionBase
{
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
  using ExpressionAdapter = ParserToGridFunctionAdapter<typename Space::EntitySet, RF>;

  //! world dimension
  static constexpr int dim = LocalBasisTraits::dimDomain;

  //! range dimension
  static constexpr int dim_range = LocalBasisTraits::dimRange;

  // this operator only support scalar ranges
  static_assert(dim_range == 1);

  ParameterTree _config;

  Space _space;

  Logging::Logger _logger;

  std::vector<std::unique_ptr<ExpressionAdapter>> _diffusion_gf;
  std::vector<std::unique_ptr<ExpressionAdapter>> _reaction_gf;
  std::vector<std::unique_ptr<ExpressionAdapter>> _jacobian_gf;
  mutable std::vector<RF> _u;

  std::vector<std::vector<std::size_t>> _component_pattern;

  /**
   * @brief      Constructs a new instance.
   *
   * @todo       Make integration order variable depending on user requirements
   *             and polynomail order of the local finite element
   *
   * @param[in]  entity_set      The entity set where this local operator is valid
   * @param[in]  config          The configuration tree
   * @param[in]  finite_element  The local finite element
   * @param[in]  rule_order      order of the quadrature rule
   */
  LocalOperatorDiffusionReactionBase(const Space& space,
                                     const ParameterTree& config)
    : _config(config)
    , _space(space)
    , _logger(Logging::Logging::componentLogger({}, "model"))
  {
    assert(_space.degree() == _config.sub("diffusion").getValueKeys().size());
    assert((_space.degree() * _space.degree()) == _config.sub("reaction.jacobian").getValueKeys().size());
    create_pattern_and_gf_expressions();
  }

  LocalOperatorDiffusionReactionBase(const LocalOperatorDiffusionReactionBase& other)
    : LocalOperatorDiffusionReactionBase(other._space, other._config)
  {}

  LocalOperatorDiffusionReactionBase(LocalOperatorDiffusionReactionBase&& other)
    : _config(std::move(other._config))
    , _space(std::move(other._space))
    , _logger(std::move(other._logger))
    , _diffusion_gf(std::move(other._diffusion_gf))
    , _reaction_gf(std::move(other._reaction_gf))
    , _jacobian_gf(std::move(other._jacobian_gf))
    , _u(std::move(other._u))
    , _component_pattern(std::move(other._component_pattern))
  {}

  LocalOperatorDiffusionReactionBase& operator=(const LocalOperatorDiffusionReactionBase& other) {
    _config = other._config;
    _space = other._space;
    _logger = other._logger;
    create_pattern_and_gf_expressions();
    return *this;
  }

  LocalOperatorDiffusionReactionBase& operator=(LocalOperatorDiffusionReactionBase&& other) {
    _config = std::move(other._config);
    _space = std::move(other._space);
    _logger = std::move(other._logger);
    _diffusion_gf = std::move(other._diffusion_gf);
    _reaction_gf = std::move(other._reaction_gf);
    _jacobian_gf = std::move(other._jacobian_gf);
    _u = std::move(other._u);
    _component_pattern = std::move(other._component_pattern);
    return *this;
  }

  void create_pattern_and_gf_expressions()
  {
    TRACE_EVENT("dune", "LocalOperator::ParserSetUp");
    _logger.trace("creating pattern and grid function expressions"_fmt);

    _u.resize(_space.degree());
    _diffusion_gf.resize(_space.degree());
    _reaction_gf.resize(_space.degree());
    _jacobian_gf.resize(_space.degree() * _space.degree());

    auto diffusion_config = _config.sub("diffusion");
    auto reaction_config = _config.sub("reaction");
    auto jacobian_config = _config.sub("reaction.jacobian");

    auto add_u = [&](auto& parser){
      for (std::size_t comp = 0; comp != _space.degree(); ++comp) {
        auto comp_space = _space.subSpace(Assembler::multiIndex(comp));
        parser.define_variable(comp_space.name(), &_u[comp]);
      }
    };

    _component_pattern.assign(_space.degree(), {});

    for (std::size_t comp_i = 0; comp_i !=  _space.degree(); ++comp_i) {
      auto comp_space_i = _space.subSpace(Assembler::multiIndex(comp_i));
      std::string var = comp_space_i.name();

      std::string d_eq = diffusion_config.template get<std::string>(var);
      std::string r_eq = reaction_config.template get<std::string>(var);

      _diffusion_gf[comp_i] = std::make_unique<ExpressionAdapter>(_space.entitySet(), make_parser());
      _diffusion_gf[comp_i]->parser().set_expression(d_eq);
      _diffusion_gf[comp_i]->parser().compile();

      _reaction_gf[comp_i] = std::make_unique<ExpressionAdapter>(_space.entitySet(), make_parser());
      _reaction_gf[comp_i]->parser().set_expression(r_eq);
      add_u(_reaction_gf[comp_i]->parser());
      _reaction_gf[comp_i]->parser().compile();

      _component_pattern[comp_i].push_back(comp_i);

      for (std::size_t comp_j = 0; comp_j != _space.degree(); ++comp_j) {
        auto comp_space_j = _space.subSpace(Assembler::multiIndex(comp_j));
        const auto jac_index = _space.degree() * comp_i + comp_j;
        std::string j_eq = jacobian_config.template get<std::string>(fmt::format("d{}__d{}", comp_space_i.name(), comp_space_j.name()));

        _jacobian_gf[jac_index] = std::make_unique<ExpressionAdapter>(_space.entitySet(), make_parser());
        _jacobian_gf[jac_index]->parser().set_expression(j_eq);
        add_u(_jacobian_gf[jac_index]->parser());
        _jacobian_gf[jac_index]->parser().compile();

        if (r_eq.find(comp_space_j.name()) != std::string::npos)
          _component_pattern[comp_i].push_back(comp_j);
      }
      std::sort(begin(_component_pattern[comp_i]), end(_component_pattern[comp_i]));
    }

    // log compartment pattern
    _logger.debug("Compartment jacobian pattern:"_fmt);
    for (std::size_t comp_i = 0; comp_i != _component_pattern.size(); ++comp_i) {
      auto comp_space_i = _space.subSpace(Assembler::multiIndex(comp_i));
      for (auto comp_j : _component_pattern[comp_i]) {
        auto comp_space_j = _space.subSpace(Assembler::multiIndex(comp_j));
        _logger.debug(2, "{} -> {}"_fmt, comp_space_i.name(), comp_space_j.name());
      }
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
        gf->setTime(t);
    for (auto& gf : _reaction_gf)
      if (gf)
        gf->setTime(t);
    for (auto& gf : _diffusion_gf)
      if (gf)
        gf->setTime(t);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_OPERATOR_DIFFUSION_REACTION_BASE_HH
