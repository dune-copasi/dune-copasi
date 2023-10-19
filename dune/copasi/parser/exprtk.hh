#ifndef DUNE_COPASI_EXPRTK_PARSER_HH
#define DUNE_COPASI_EXPRTK_PARSER_HH

// #include <dune/logging/logging.hh>

#include <dune/copasi/parser/parser.hh>

#include <exprtk.hpp>

#include <string>
#include <vector>
#include <unordered_map>

namespace Dune::Copasi {


class ExprTkParser final : public Parser {
  static const std::size_t max_functions = 100;

  template<class F>
  struct FunctionID {
    std::shared_ptr<const exprtk::parser<RangeField>> parser;
    std::string name;
    F function;
  };

public:
  ExprTkParser()
    : Parser{}
    , _parser{std::make_shared<exprtk::parser<RangeField>>()}
  {}

  ~ExprTkParser() override {};

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
    _symbol_table.add_constant(symbol, value);
  };

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
    return ExprTkParser::_function2d[i]->function(arg0,arg1);
  }

  template<std::size_t i>
  static RangeField function_wrapper_3d(RangeField arg0, RangeField arg1, RangeField arg2)
  {
    assert(ExprTkParser::_function3d[i]);
    return ExprTkParser::_function3d[i]->function(arg0,arg1,arg2);
  }

  void compile() override {
    Parser::compile();
    assert(not _compiled);
    assert(size(_symbols) == size(_variables));

    for (std::size_t i = 0; i < size(_symbols); ++i) {
      _symbol_table.add_variable(_symbols[i],*_variables[i]);
    }

    auto indices = Dune::range(std::integral_constant<std::size_t, max_functions>{});
    Dune::Hybrid::forEach(indices, [&](auto i) {
      auto i0 = _function0d.find(i);
      auto i1 = _function1d.find(i);
      auto i2 = _function2d.find(i);
      auto i3 = _function3d.find(i);

      if (i0 != end(_function0d) and (i0->second) and (i0->second->parser == _parser))
        _symbol_table.add_function(i0->second->name, function_wrapper_0d<i>);
      if (i1 != end(_function1d) and (i1->second) and (i1->second->parser == _parser))
        _symbol_table.add_function(i1->second->name, function_wrapper_1d<i>);
      if (i2 != end(_function2d) and (i2->second) and i2->second->parser == _parser)
        _symbol_table.add_function(i2->second->name, function_wrapper_2d<i>);
      if (i3 != end(_function3d) and (i3->second) and i3->second->parser == _parser)
        _symbol_table.add_function(i3->second->name, function_wrapper_3d<i>);
    });

    _tk_expression.register_symbol_table(_symbol_table);
    _parser->compile(_expression,_tk_expression);
    _compiled = true;
  }

  [[nodiscard]] RangeField eval() const noexcept override {
    assert(_compiled);
    return _tk_expression.value();
  }

private:
  using Parser::_compiled;
  using Parser::_expression;
  using Parser::_symbols;
  using Parser::_variables;

  std::shared_ptr<exprtk::parser<RangeField>> _parser;
  exprtk::expression<RangeField> _tk_expression;
  exprtk::symbol_table<RangeField> _symbol_table;

  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID0D>> _function0d;
  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID1D>> _function1d;
  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID2D>> _function2d;
  static std::unordered_map<std::size_t,std::unique_ptr<FunctionID3D>> _function3d;
};

std::unordered_map<std::size_t,std::unique_ptr<ExprTkParser::FunctionID0D>> ExprTkParser::_function0d = {};
std::unordered_map<std::size_t,std::unique_ptr<ExprTkParser::FunctionID1D>> ExprTkParser::_function1d = {};
std::unordered_map<std::size_t,std::unique_ptr<ExprTkParser::FunctionID2D>> ExprTkParser::_function2d = {};
std::unordered_map<std::size_t,std::unique_ptr<ExprTkParser::FunctionID3D>> ExprTkParser::_function3d = {};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_EXPRTK_PARSER_HH
