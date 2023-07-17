#ifndef DUNE_COPASI_EXPRTK_PARSER_HH
#define DUNE_COPASI_EXPRTK_PARSER_HH

#include <dune/copasi/parser/parser.hh>

#include <exprtk.hpp>

#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace Dune::Copasi {

class ExprTkParser final : public Parser
{
  static const std::size_t max_functions = 100;

  template<class F>
  struct FunctionID
  {
    std::shared_ptr<const exprtk::parser<RangeField>> parser;
    std::string name;
    F function;
  };

public:
  ExprTkParser();

  ~ExprTkParser() override final;

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
  void define_constant(const std::string& symbol, const RangeField& value) override final;

  void define_function(const std::string& symbol, const Function0D& function) override final;

  void define_function(const std::string& symbol, const Function1D& function) override final;

  void define_function(const std::string& symbol, const Function2D& function) override final;

  void define_function(const std::string& symbol, const Function3D& function) override final;

  // private: ????

  template<std::size_t i>
  static RangeField function_wrapper_0d()
  {
    assert(ExprTkParser::_function0d[i]);
    return ExprTkParser::_function0d[i]->function();
  }

  template<std::size_t i>
  static RangeField function_wrapper_1d(RangeField arg0)
  {
    assert(ExprTkParser::_function1d[i]);
    return ExprTkParser::_function1d[i]->function(arg0);
  }

  template<std::size_t i>
  static RangeField function_wrapper_2d(RangeField arg0, RangeField arg1)
  {
    assert(ExprTkParser::_function2d[i]);
    return ExprTkParser::_function2d[i]->function(arg0, arg1);
  }

  template<std::size_t i>
  static RangeField function_wrapper_3d(RangeField arg0, RangeField arg1, RangeField arg2)
  {
    assert(ExprTkParser::_function3d[i]);
    return ExprTkParser::_function3d[i]->function(arg0, arg1, arg2);
  }

  void compile() override final;

  [[nodiscard]] RangeField operator()() const noexcept override final;

private:
  void register_functions();
  void unregister_functions();

  using Parser::_compiled;
  using Parser::_expression;
  using Parser::_symbols;
  using Parser::_variables;

  std::shared_ptr<exprtk::parser<RangeField>> _parser;
  exprtk::expression<RangeField> _tk_expression;
  exprtk::symbol_table<RangeField> _symbol_table;

  static inline std::map<std::size_t, std::unique_ptr<FunctionID0D>> _function0d = {};
  static inline std::map<std::size_t, std::unique_ptr<FunctionID1D>> _function1d = {};
  static inline std::map<std::size_t, std::unique_ptr<FunctionID2D>> _function2d = {};
  static inline std::map<std::size_t, std::unique_ptr<FunctionID3D>> _function3d = {};

  // recursive mutex are needed because functions may hold parsers that need to be destructed
  static inline std::recursive_mutex _mutex = {};
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_EXPRTK_PARSER_HH
