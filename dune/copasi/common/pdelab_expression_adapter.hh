// #ifndef DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
// #define DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH

// #include <dune/pdelab/common/function.hh>

// #include <dune/logging/logging.hh>

// #include <dune/common/dynvector.hh>
// #include <dune/common/parametertree.hh>

// #include <muParser.h>

// #include <string>
// #include <type_traits>
// #include <memory>

// namespace Dune::Copasi {

// namespace Impl {
//   /**
//    * @brief      Output information on the parser error and throw DUNE exception
//    *
//    * @param[in]  e     Exception thrown by the parser
//    */
//   template<class E>
//   void handle_parser_error(const E& e)
//   {
//     DUNE_THROW(IOError,
//     "Evaluating muParser expression failed:" << std::endl
//     << "  Parsed expression:   " << e.GetExpr() << std::endl
//     << "  Token:               " << e.GetToken() << std::endl
//     << "  Error position:      " << e.GetPos() << std::endl
//     << "  Error code:          " << int(e.GetPos()) << std::endl
//     << "  Error message:       " << e.GetMsg() << std::endl
//     );
//   }
// }

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
//                                   const std::vector<std::string>& other_variables = {})
//     : _logger(Logging::Logging::componentLogger({}, "model"))
//     , _gv(grid_view)
//     , _time(0.)
//     , _other_value(other_variables.size())
//     , _expr(equation)
//     , _compiled(false)
//   {

//     constexpr int dim = Traits::dimDomain;

//     _logger.trace("initialize parser with constant variables"_fmt);
//     _parser.DefineConst("pi", StandardMathematicalConstants<double>::pi());
//     _parser.DefineConst("dim", dim);

//     _parser.DefineVar("t", &_time);
//     _parser.DefineVar("x", &_pos_global[0]);
//     _parser.DefineVar("y", &_pos_global[1]);

//     if constexpr (dim == 3)
//       _parser.DefineVar("z", &_pos_global[2]);

//     for (size_t i = 0; i < other_variables.size(); i++) {
//       _logger.trace("define extra variable: {}"_fmt, other_variables[i]);
//       _parser.DefineVar(other_variables[i], &_other_value[i]);
//     }

//     if (do_compile_parser)
//       compile_parser();
//   }

//   ExpressionToGridFunctionAdapter(const ExpressionToGridFunctionAdapter& other)
//     : _logger(Logging::Logging::componentLogger({}, "model"))
//     , _gv(other._gv)
//     , _time(other._time)
//     , _other_value(other._other_value)
//     , _parser(other._parser)
//     , _expr(other._expr)
//     , _compiled(false)
//   {
//     if (other._compiled)
//       compile_parser();
//   }

//   ExpressionToGridFunctionAdapter(ExpressionToGridFunctionAdapter&& other)
//     : _logger(Logging::Logging::componentLogger({}, "model"))
//     , _gv(other._gv)
//     , _time(other._time)
//     , _other_value(std::move(other._other_value))
//     , _parser(std::move(other._parser))
//     , _expr(std::move(other._expr))
//     , _compiled(false)
//   {
//     if (other._compiled)
//       compile_parser();
//   }

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
//   void set_time(double t) { _time = t; }

// private:

//   Logging::Logger _logger;

//   GV _gv;

//   /// Cache for the evaluation position
//   mutable typename Traits::DomainType _pos_global;

//   double _time;

//   DynamicVector<RF> _other_value;

//   /// The parser instance
//   mu::Parser _parser;

//   std::string _expr;
//   bool _compiled;
// };

// /**
//  * @brief      Gets the muparser expressions.
//  *
//  * @param[in]  expressions_config  The expressions configuration
//  * @param[in]  gf_grid_view        The grid view for the grid function
//  * @param[in]  compile             True to compile expression at construction
//  *
//  * @tparam     GFGridView          Grid view type
//  * @tparam     RF                  Range field type
//  *
//  * @return     Vector with muparser expressions pointers
//  */
// template<class GFGridView, class RF = double>
// auto get_muparser_expressions(
//   const ParameterTree& expressions_config,
//   const GFGridView& gf_grid_view,
//   bool compile = true)
// {
//   const auto& vars = expressions_config.getValueKeys();

//   using GridFunction = ExpressionToGridFunctionAdapter<GFGridView, RF>;
//   std::vector<std::shared_ptr<GridFunction>> functions;

//   for (std::size_t i = 0; i < vars.size(); i++) {
//     auto gf = std::make_shared<GridFunction>(
//       gf_grid_view, expressions_config[vars[i]], compile);
//     functions.emplace_back(gf);
//     assert(functions.size() == i + 1);
//     if (compile)
//       functions[i]->compile_parser();
//   }

//   return functions;
// }

// } // namespace Dune::Copasi

// #endif // DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
