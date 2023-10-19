#ifndef DUNE_COPASI_MU_PARSER_HH
#define DUNE_COPASI_MU_PARSER_HH

#include <dune/copasi/parser/parser.hh>

// #include <dune/logging/logging.hh>

#include <dune/common/hybridutilities.hh>
#include <dune/common/exceptions.hh>

#include <muParser.h>

#include <string>
#include <memory>

namespace Dune::Copasi {


class MuParser final : public Parser {
  static const std::size_t max_functions = 100;

  template<class F>
  struct FunctionID {
    std::shared_ptr<const mu::Parser> parser;
    std::string name;
    F function;
  };

public:
  MuParser()
    : Parser{}
    , _parser{std::make_shared<mu::Parser>()}
  {
    _parser->DefineConst("pi", StandardMathematicalConstants<double>::pi());
  }

  ~MuParser() override {};

  using RangeField = typename Parser::RangeField;

  using Function0D = typename Parser::Function0D;
  using Function1D = typename Parser::Function1D;
  using Function2D = typename Parser::Function2D;
  using Function3D = typename Parser::Function3D;

private:
  using FunctionID0D = FunctionID<Function0D>;
  using FunctionID1D = FunctionID<Function1D>;
  using FunctionID2D = FunctionID<Function2D>;
  using FunctionID3D = FunctionID<Function3D>;

public:

  void define_constant(const std::string& symbol, const RangeField& value) override {
    _parser->DefineConst(symbol, value);
  }

  void define_function(const std::string& symbol, Function0D&& function) override {
    for (std::size_t i = 0; i < max_functions; ++i)
      if (not _function0d[i]) {
        _function0d[i] = std::make_unique<FunctionID0D>(FunctionID0D{_parser, symbol, std::move(function) });
        return;
      }
    DUNE_THROW(Exception,""); // not enough instantiated functions
  }

  void define_function(const std::string& symbol, Function1D&& function) override {
    for (std::size_t i = 0; i < max_functions; ++i)
      if (not _function1d[i]) {
        _function1d[i] = std::make_unique<FunctionID1D>(FunctionID1D{_parser, symbol, std::move(function) });
        return;
      }
    DUNE_THROW(Exception,""); // not enough instantiated functions
  }

  void define_function(const std::string& symbol, Function2D&& function) override {
    for (std::size_t i = 0; i < max_functions; ++i)
      if (not _function2d[i]) {
        _function2d[i] = std::make_unique<FunctionID2D>(FunctionID2D{_parser, symbol, std::move(function) });
        return;
      }
    DUNE_THROW(Exception,""); // not enough instantiated functions
  }

  void define_function(const std::string& symbol, Function3D&& function) override {
    for (std::size_t i = 0; i < max_functions; ++i)
      if (not _function3d[i]) {
        _function3d[i] = std::make_unique<FunctionID3D>(FunctionID3D{_parser, symbol, std::move(function) });
        return;
      }
    DUNE_THROW(Exception,""); // not enough instantiated functions
  }

// private: ????

  template<std::size_t i>
  static RangeField function_wrapper_0d()
  {
    assert(MuParser::_function0d[i]);
    return MuParser::_function0d[i]->function();
  }

  template<std::size_t i>
  static RangeField function_wrapper_1d(RangeField arg0)
  {
    assert(MuParser::_function1d[i]);
    return MuParser::_function1d[i]->function(arg0);
  }

  template<std::size_t i>
  static RangeField function_wrapper_2d(RangeField arg0, RangeField arg1)
  {
    assert(MuParser::_function2d[i]);
    return MuParser::_function2d[i]->function(arg0,arg1);
  }

  template<std::size_t i>
  static RangeField function_wrapper_3d(RangeField arg0, RangeField arg1, RangeField arg2)
  {
    assert(MuParser::_function3d[i]);
    return MuParser::_function3d[i]->function(arg0,arg1,arg2);
  }

public:
  void compile() override {
    Parser::compile();
    assert(not _compiled);
    assert(size(_symbols) == size(_variables));

    for (std::size_t i = 0; i < size(_symbols); ++i)
      _parser->DefineVar(_symbols[i], _variables[i]);

    auto indices = Dune::range(std::integral_constant<std::size_t, max_functions>{});
    Dune::Hybrid::forEach(indices, [&](auto i) {
      auto i0 = _function0d.find(i);
      auto i1 = _function1d.find(i);
      auto i2 = _function2d.find(i);
      auto i3 = _function3d.find(i);

      if (i0 != end(_function0d) and (i0->second) and (i0->second->parser == _parser))
        _parser->DefineFun(i0->second->name, function_wrapper_0d<i>);
      if (i1 != end(_function1d) and (i1->second) and (i1->second->parser == _parser))
        _parser->DefineFun(i1->second->name, function_wrapper_1d<i>);
      if (i2 != end(_function2d) and (i2->second) and i2->second->parser == _parser)
        _parser->DefineFun(i2->second->name, function_wrapper_2d<i>);
      if (i3 != end(_function3d) and (i3->second) and i3->second->parser == _parser)
        _parser->DefineFun(i3->second->name, function_wrapper_3d<i>);
    });
    _parser->SetExpr(_expression);
    _compiled = true;

    try {
      _parser->Eval();
    } catch (mu::Parser::exception_type& e) {
      DUNE_THROW(IOError,
      "Evaluating muParser expression failed:" << std::endl
      << "  Parsed expression:   " << e.GetExpr() << std::endl
      << "  Token:               " << e.GetToken() << std::endl
      << "  Error position:      " << e.GetPos() << std::endl
      << "  Error code:          " << int(e.GetPos()) << std::endl
      << "  Error message:       " << e.GetMsg() << std::endl
      );
    }
  }

  [[nodiscard]] RangeField eval() const noexcept override {
    try {
      assert(_compiled);
      return _parser->Eval();
    } catch (mu::Parser::exception_type& e) {
      std::cerr
        << "Evaluating muParser expression failed:" << std::endl
        << "  Parsed expression:   " << e.GetExpr() << std::endl
        << "  Token:               " << e.GetToken() << std::endl
        << "  Error position:      " << e.GetPos() << std::endl
        << "  Error code:          " << int(e.GetPos()) << std::endl
        << "  Error message:       " << e.GetMsg() << std::endl;
      std::terminate();
    } catch (...) {
      std::terminate();
    }
  }

private:
  using Parser::_compiled;
  using Parser::_expression;
  using Parser::_symbols;
  using Parser::_variables;

  std::shared_ptr<mu::Parser> _parser;

  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID0D>> _function0d;
  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID1D>> _function1d;
  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID2D>> _function2d;
  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID3D>> _function3d;
};

std::unordered_map<std::size_t,std::unique_ptr<MuParser::FunctionID0D>> MuParser::_function0d = {};
std::unordered_map<std::size_t,std::unique_ptr<MuParser::FunctionID1D>> MuParser::_function1d = {};
std::unordered_map<std::size_t,std::unique_ptr<MuParser::FunctionID2D>> MuParser::_function2d = {};
std::unordered_map<std::size_t,std::unique_ptr<MuParser::FunctionID3D>> MuParser::_function3d = {};




// /**
//  * @brief      Converts an interface to match an expression to a PDELab grid
//  *             function.
//  * @details    The resulting grid view is only for scalar expressions
//  *
//  * @tparam     GV    Grid View
//  * @tparam     RF    Range Field
//  */
// template<typename GV, typename RF>
// class ExpressionToGridFunctionAdapter
//   : public PDELab::GridFunctionBase<
//       PDELab::GridFunctionTraits<GV, RF, 1, FieldVector<RF, 1>>,
//       ExpressionToGridFunctionAdapter<GV, RF>>
// {
//   void configure(const bool& do_compile_parser){
//     constexpr int dim = Traits::dimDomain;

//     _logger.trace("initialize parser with constant variables"_fmt);
//     _parser.DefineConst("pi", StandardMathematicalConstants<double>::pi());
//     _parser.DefineConst("dim", dim);

//     _parser.DefineVar("t", &_time);
//     _parser.DefineVar("x", &_pos_global[0]);
//     _parser.DefineVar("y", &_pos_global[1]);

//     if constexpr (dim == 3)
//       _parser.DefineVar("z", &_pos_global[2]);

//     for (size_t i = 0; i < _other_variables_name.size(); i++) {
//       _logger.trace("define extra variable: {}"_fmt, _other_variables_name[i]);
//       _parser.DefineVar(_other_variables_name[i], &_other_value[i]);
//     }

//     if (do_compile_parser)
//       compile_parser();
//   }

// public:
//   using Traits = PDELab::GridFunctionTraits<GV, RF, 1, FieldVector<RF, 1>>;

//   /**
//    * @brief      Constructs a new instance.
//    *
//    * @param[in]  grid_view          The grid view
//    * @param[in]  equation           The math expression
//    * @param[in]  do_compile_parser  Bool to compile parser at object
//    * construction
//    * @param[in]  other_variables    Extra varialbes names to be available in the
//    * expression
//    */
//   ExpressionToGridFunctionAdapter(const GV& grid_view,
//                                   const std::string& equation,
//                                   bool do_compile_parser = true,
//                                   const std::vector<std::string>& other_variables_name = {})
//     : _logger(Logging::Logging::componentLogger({}, "model"))
//     , _gv(grid_view)
//     , _time(0.)
//     , _other_variables_name(other_variables_name)
//     , _other_value(_other_variables_name.size())
//     , _expr(equation)
//     , _compiled(false)
//   {
//     configure(do_compile_parser);
//   }

//   ExpressionToGridFunctionAdapter(const ExpressionToGridFunctionAdapter& other)
//     : _logger(other._logger)
//     , _gv(other._gv)
//     , _time(other._time)
//     , _other_value(other._other_value)
//     , _other_variables_name(other._other_variables_name)
//     , _expr(other._expr)
//     , _compiled(false)
//   {
//     configure(other._compiled);
//   }

// public:
//   /**
//    * @brief      Gets a reference to the grid view
//    *
//    * @return     The grid view.
//    */
//   inline const GV& getGridView() const { return _gv; }

//   /**
//    * @brief      Evaluates extended function on a element
//    *
//    * @param[in]  e     Entity to operate with
//    * @param[in]  x     Local coordinates in the entity
//    * @param      y     Resulting value
//    *
//    * @tparam     E     Entity
//    * @tparam     D     Domain
//    * @tparam     R     Range
//    */
//   template<class E, class D, class R>
//   void evaluate(const E& e, const D& x, R& y) const
//   {
//     assert(_compiled);
//     // update position storage
//     _pos_global = e.geometry().global(x);

//     // evaluate the expression (with new position)
//     try {
//       y = _parser.Eval();
//     } catch (mu::Parser::exception_type& e) {
//       Impl::handle_parser_error(e);
//     }
//   }

//   /**
//    * @brief      Updates the given other values set at construction.
//    * @details    The extra values have to have the same order as entred in the
//    *             object construction
//    *
//    * @param[in]  other_value  The other values
//    */
//   inline void update(const DynamicVector<RF>& other_value)
//   {
//     assert(_other_value.size() == other_value.size());
//     _other_value = other_value;
//   }

//   /**
//    * @brief      Compiles the parser and checks that it is able to be evaluated
//    */
//   void compile_parser()
//   {
//     try {
//       _logger.trace("compile expression: {}"_fmt, _expr);
//       _parser.SetExpr(_expr);
//       // try to evaluate once
//       _parser.Eval();
//     } catch (mu::Parser::exception_type& e) {
//       Impl::handle_parser_error(e);
//     }
//     _compiled = true;
//   }

//   /**
//    * @brief      Get parser
//    *
//    * @return     Reference to internal parser
//    */
//   mu::Parser& parser()
//   {
//     assert(not _compiled);
//     return _parser;
//   }

//   /**
//    * @brief      Sets the time.
//    *
//    * @param[in]  t     The new time
//    */
//   void setTime(double t) { _time = t; }

// private:

//   Logging::Logger _logger;

//   GV _gv;

//   /// Cache for the evaluation position
//   mutable typename Traits::DomainType _pos_global;

//   double _time;

//   std::vector<std::string> _other_variables_name;
//   DynamicVector<RF> _other_value;

//   /// The parser instance
//   mu::Parser _parser;

//   std::string _expr;
//   bool _compiled;
// };

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MU_PARSER_HH
