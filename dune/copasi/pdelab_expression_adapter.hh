#ifndef DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
#define DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH

#include <dune/copasi/tiff_grayscale.hh>

#include <dune/pdelab/common/function.hh>

#include <dune/logging/logging.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/parametertree.hh>

#include <muParser.h>

#include <algorithm>
#include <string>
#include <type_traits>

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
                                  bool do_compile_parser = true,
                                  std::vector<std::string> other_variables = {})
    : _logger(Logging::Logging::componentLogger({}, "model"))
    , _gv(grid_view)
    , _time(0.)
    , _other_value(other_variables.size())
    , _expr(equation)
    , _compiled(false)
  {

    constexpr int dim = Traits::dimDomain;

    std::sort(other_variables.begin(), other_variables.end());

    _logger.trace("initialize parser with constant variables"_fmt);
    _parser.DefineConst("pi", StandardMathematicalConstants<double>::pi());
    _parser.DefineConst("dim", dim);

    _parser.DefineVar("t", &_time);
    _parser.DefineVar("x", &_pos_global[0]);
    _parser.DefineVar("y", &_pos_global[1]);

    if constexpr (dim == 3)
      _parser.DefineVar("z", &_pos_global[2]);

    for (size_t i = 0; i < other_variables.size(); i++) {
      _logger.trace("define extra variable: {}"_fmt, other_variables[i]);
      _parser.DefineVar(other_variables[i], &_other_value[i]);
    }

    if (do_compile_parser)
      compile_parser();

    _logger.debug("ExpressionToGridFunctionAdapter constructed"_fmt);
  }

  ~ExpressionToGridFunctionAdapter()
  {
    _logger.debug("ExpressionToGridFunctionAdapter deconstructed"_fmt);
  }

  //! get a reference to the grid view
  inline const GV& getGridView() const { return _gv; }

  //! evaluate extended function on element
  template<class E, class D, class R>
  void evaluate(const E& e, const D& x, R& y) const
  {
    assert(_compiled);
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

  void compile_parser()
  {
    try {
      _logger.trace("compile expression: {}"_fmt, _expr);
      _parser.SetExpr(_expr);
      // try to evaluate once
      _parser.Eval();
    } catch (mu::Parser::exception_type& e) {
      handle_parser_error(e);
    }
    _compiled = true;
  }

  mu::Parser& parser()
  {
    assert(not _compiled);
    return _parser;
  }

  void set_time(double t) { _time = t; }

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

  double _time;

  DynamicVector<RF> _other_value;

  /// The parser instance
  mu::Parser _parser;

  std::string _expr;
  bool _compiled;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH