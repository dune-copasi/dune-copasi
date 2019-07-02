#ifndef DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
#define DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH

#include <dune/pdelab/common/function.hh>

#include <dune/logging/logging.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/parametertree.hh>

#include <muParser.h>

#include <string>

namespace Dune::Copasi {

template<typename GV, typename RF>
class ExpressionToGridFunctionAdapter 
  : public PDELab::GridFunctionBase<PDELab::GridFunctionTraits<GV,RF,-1,DynamicVector<RF>>,
                                    ExpressionToGridFunctionAdapter<GV,RF> >
{
public:
  using Traits = PDELab::GridFunctionTraits<GV,RF,-1,DynamicVector<RF>>;

  //! construct from grid view
  ExpressionToGridFunctionAdapter (const GV& grid_view, const ParameterTree& config, const ParameterTree& extra_config = {}) 
    : _logger(Logging::Logging::componentLogger(config,"default"))
    , _gv(grid_view)
    , _size(config.getValueKeys().size())
    , _extra_var(extra_config.getValueKeys().size())
    , _parser(_size)
  {
    assert(_size > 0);

    constexpr int dim = Traits::dimDomain;

    const auto& keys = config.getValueKeys();

    for (int i = 0; i < keys.size(); ++i)
    {
      _logger.trace("setting up variable: {}"_fmt, keys[i]);

      _logger.trace("initialize parser with constant variables"_fmt);
      _parser[i].DefineConst("pi",StandardMathematicalConstants<double>::pi());
      _parser[i].DefineConst("dim", dim);

      _parser[i].DefineVar("x", &_pos_global[0]);
      _parser[i].DefineVar("y", &_pos_global[1]);
      if constexpr (dim == 3)
        _parser[i].DefineVar("z", &_pos_global[2]);

      const auto& extra_keys = extra_config.getValueKeys();
      for (int j = 0; j < extra_keys.size(); ++j)
      {
        _logger.trace("define extra variable: {}"_fmt, extra_keys[j]);
        _parser[i].DefineVar(extra_keys[j], &_extra_var[j]);
      }


      // set up parser expression
      try {
        std::string equation = config.template get<std::string>(keys[i]);
        _logger.trace("compile expression: {}"_fmt, equation);
        _parser[i].SetExpr(equation);
        // try to evaluate once
        _parser[i].Eval();
      } catch (mu::Parser::exception_type& e) {
        handle_parser_error(e);
      }
    }

    _logger.debug("ExpressionToGridFunctionAdapter constructed"_fmt);
  }

  ~ExpressionToGridFunctionAdapter()
  {
    _logger.debug("ExpressionToGridFunctionAdapter destructed"_fmt);
  }

  //! get a reference to the grid view
  inline const GV& getGridView () const {return _gv;}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y.resize(_size);
    // update position storage
    _pos_global = e.geometry().global(x);

    // evaluate the expression (with new position)
    try {
      for (int i = 0; i < _size; ++i)
        y[i] = _parser[i].Eval();
    }
    catch (mu::Parser::exception_type& e) {
      handle_parser_error(e);
    }
  }

  void bind(const typename Traits::ElementType& e, const DynamicVector<RF>& extra_var)
  {
    _extra_var = extra_var;
  }

private:
  /// Output information on the parser error and throw DUNE exception
  /**
   *  \param e Exception thrown by the parser
   *  \throw IOError (always throws)
   */
  void handle_parser_error (const mu::Parser::exception_type& e) const
  {
    _logger.error("Evaluating analytic initial condition failed:"_fmt);
    _logger.error("  Parsed expression:   {}"_fmt,e.GetExpr());
    _logger.error("  Token:               {}"_fmt,e.GetToken());
    _logger.error("  Error position:      {}"_fmt,e.GetPos());
    _logger.error("  Error code:          {}"_fmt,int(e.GetCode()));
    _logger.error("  Error message:       {}"_fmt,e.GetMsg());
    DUNE_THROW(IOError, "Error evaluating analytic initial condition");
  }

  Logging::Logger  _logger;

  GV _gv;
  const std::size_t _size;

  /// Cache for the evaluation position
  mutable typename Traits::DomainType _pos_global;

  DynamicVector<RF> _extra_var;

  /// The parser instance
  std::vector<mu::Parser> _parser;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH