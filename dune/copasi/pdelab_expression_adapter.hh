#ifndef DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
#define DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH

#include <dune/pdelab/common/function.hh>

#include <dune/logging/logging.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/parametertree.hh>

#include <muParser.h>

#include <algorithm>
#include <string>

namespace Dune::Copasi {

template<typename GV, typename RF>
class ExpressionToGridFunctionAdapter
  : public PDELab::GridFunctionBase<
      PDELab::GridFunctionTraits<GV, RF, 1, FieldVector<RF, 1>>,
      ExpressionToGridFunctionAdapter<GV, RF>>
{
public:
  using Traits = PDELab::GridFunctionTraits<GV, RF, 1, FieldVector<RF, 1>>;

  //! construct from grid view
  ExpressionToGridFunctionAdapter(const GV& grid_view,
                                  const std::string& equation,
                                  std::vector<std::string> other_variables = {})
    : _logger(Logging::Logging::componentLogger({}, "default"))
    , _gv(grid_view)
    , _other_value(other_variables.size())
  {

    constexpr int dim = Traits::dimDomain;

    std::sort(other_variables.begin(), other_variables.end());

    _logger.trace("initialize parser with constant variables"_fmt);
    _parser.DefineConst("pi", StandardMathematicalConstants<double>::pi());
    _parser.DefineConst("dim", dim);

    _parser.DefineVar("x", &_pos_global[0]);
    _parser.DefineVar("y", &_pos_global[1]);
    if constexpr (dim == 3)
      _parser.DefineVar("z", &_pos_global[2]);

    // set up parser expression
    try {
      for (size_t i = 0; i < other_variables.size(); i++) {
        _logger.trace("define extra variable: {}"_fmt, other_variables[i]);
        _parser.DefineVar(other_variables[i], &_other_value[i]);
      }
      _logger.trace("compile expression: {}"_fmt, equation);
      _parser.SetExpr(equation);
      // try to evaluate once
      _parser.Eval();
    } catch (mu::Parser::exception_type& e) {
      handle_parser_error(e);
    }

    _logger.debug("ExpressionToGridFunctionAdapter constructed"_fmt);
  }

  ~ExpressionToGridFunctionAdapter()
  {
    _logger.debug("ExpressionToGridFunctionAdapter destructed"_fmt);
  }

  //! get a reference to the grid view
  inline const GV& getGridView() const { return _gv; }

  //! evaluate extended function on element
  template<class E, class D, class R>
  void evaluate(const E& e, const D& x, R& y) const
  {
    // update position storage
    _pos_global = e.geometry().global(x);

    // evaluate the expression (with new position)
    try {
      y = _parser.Eval();
    } catch (mu::Parser::exception_type& e) {
      handle_parser_error(e);
    }
  }

  inline void update(const DynamicVector<RF>& other_value)
  {
    assert(_other_value.size() == other_value.size());
    _other_value = other_value;
  }

private:
  /// Output information on the parser error and throw DUNE exception
  /**
   *  \param e Exception thrown by the parser
   *  \throw IOError (always throws)
   */
  void handle_parser_error(const mu::Parser::exception_type& e) const
  {
    _logger.error("Evaluating analytic initial condition failed:"_fmt);
    _logger.error("  Parsed expression:   {}"_fmt, e.GetExpr());
    _logger.error("  Token:               {}"_fmt, e.GetToken());
    _logger.error("  Error position:      {}"_fmt, e.GetPos());
    _logger.error("  Error code:          {}"_fmt, int(e.GetCode()));
    _logger.error("  Error message:       {}"_fmt, e.GetMsg());
    DUNE_THROW(IOError, "Error evaluating analytic initial condition");
  }

  Logging::Logger _logger;

  GV _gv;

  /// Cache for the evaluation position
  mutable typename Traits::DomainType _pos_global;

  DynamicVector<RF> _other_value;

  /// The parser instance
  mu::Parser _parser;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH