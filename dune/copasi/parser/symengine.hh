#ifndef DUNE_COPASI_SYMENGINE_PARSER_HH
#define DUNE_COPASI_SYMENGINE_PARSER_HH

// #include <dune/logging/logging.hh>

#include <dune/copasi/parser/parser.hh>

#include <symengine/lambda_double.h>
#include <symengine/parser.h>
#include <symengine/parser/sbml/sbml_parser.h>

#include <string>
#include <functional>
#include <unordered_map>
#include <vector>

namespace Dune::Copasi {

class SymEngineParser final : public Parser
{

public:
  enum class Type {Native, SBML};

  using RangeField = typename Parser::RangeField;

  using Function0D = typename Parser::Function0D;
  using Function1D = typename Parser::Function1D;
  using Function2D = typename Parser::Function2D;
  using Function3D = typename Parser::Function3D;

  SymEngineParser(const std::string& expression, Type parser_type = Type::Native)
    : Parser{ expression }
  {
    if (Type::Native == parser_type)
      _se_expression = SymEngine::parse(_expression);
    else
      _se_expression = SymEngine::parse_sbml(_expression);
  }

  void define_constant(const std::string& symbol, const RangeField& value) override {
    _input.push_back(value);
    _const_symbols.push_back(symbol);
  }

  auto setup_function_symbol(const std::string& symbol){
      for (const auto &basic_func : function_symbols(*_se_expression)) {
        auto func = SymEngine::rcp_dynamic_cast<const SymEngine::FunctionSymbol>(basic_func);
        auto fname = func->get_name();
        if (fname == symbol) {
          SymEngine::map_basic_basic replace_map;
          replace_map[func] = SymEngine::symbol("__" + fname);
          _se_expression = _se_expression->xreplace(replace_map);
          define_constant("__" + fname, 0.);
          // we may need to delay the argument index deduction until setup time
          std::vector<std::shared_ptr<std::size_t>> args_i;
          args_i.emplace_back(std::make_shared<std::size_t>(size(_input) - 1));
          const auto& args = func->get_args();
          std::size_t count = 0;
          for (const auto &arg : args) {
            auto& arg_i = args_i.emplace_back(std::make_shared<std::size_t>(9999999));
            if (has_symbol(*_se_expression,*arg)) {
              _setup.emplace_back([this, arg_i, arg]{
                auto arg_name = SymEngine::rcp_dynamic_cast<const SymEngine::Symbol>(arg)->get_name();
                auto it = find(begin(_symbols), end(_symbols), arg_name);
                if (it != end(_symbols)){
                  (*arg_i) = distance(begin(_symbols), it);
                } else {
                  // I am not sure what has_symbol exactly does, so we may get here??
                  DUNE_THROW(Exception, "Symbol not found in main expression: " << arg_name);
                }
              });
              continue;
            } else {
              auto arg_name = "__" + fname + "__arg" + std::to_string(count++);
              define_constant(arg_name, 0.);
              *arg_i = size(_input) - 1;
              replace_map[arg] = SymEngine::symbol(arg_name);
              auto& visitor = _visitors.emplace_back();
              _setup.emplace_back([this, &visitor, arg]{
                visitor.init(_arguments, *arg);
              });
              _callbacks.emplace_back([&visitor, arg_i, this](){
                // actual update of argument value
                _input[*arg_i] = visitor.call(_input);
              });
            }
          }
          return args_i;
        }
      }
      DUNE_THROW(IOError, "Function symbol is not contained in expression!");
  }

  void define_function(const std::string& symbol,
                       Function0D&& function) override {
      auto args_i = setup_function_symbol(symbol);

      if (size(args_i) != 1)
        DUNE_THROW(IOError, "Function arguments do not match with defined function");

      _callbacks.emplace_back([args_i, this, function = std::move(function)](){
        // update function return value
        _input[*args_i[0]] = function();
      });
  };

  void define_function(const std::string& symbol,
                       Function1D&& function) override {
      auto args_i = setup_function_symbol(symbol);

      if (size(args_i) != 2)
        DUNE_THROW(IOError, "Function arguments do not match with defined function");

      _callbacks.emplace_back([args_i, this, function = std::move(function)](){
        // update function return value
        _input[*args_i[0]] = function(_input[*args_i[1]]);
      });
  };

  void define_function(const std::string& symbol,
                       Function2D&& function) override{
      auto args_i = setup_function_symbol(symbol);

      if (size(args_i) != 3)
        DUNE_THROW(IOError, "Function arguments do not match with defined function");

      _callbacks.emplace_back([args_i, this, function = std::move(function)](){
        // update function return value
        _input[*args_i[0]] = function(_input[*args_i[1]], _input[*args_i[2]]);
      });
  };

  void define_function(const std::string& symbol,
                       Function3D&& function) override{
      auto args_i = setup_function_symbol(symbol);

      if (size(args_i) != 4)
        DUNE_THROW(IOError, "Function arguments do not match with defined function");

      _callbacks.emplace_back([args_i, this, function = std::move(function)](){
        // update function return value
        _input[*args_i[0]] = function(_input[*args_i[1]], _input[*args_i[2]], _input[*args_i[3]]);
      });
  };

  void compile() override
  {
    Parser::compile();
    assert(not _compiled);
    assert(size(_symbols) == size(_variables));

    try {
      _symbols.insert(begin(_symbols),begin(_const_symbols),end(_const_symbols));
      _arguments.resize(size(_symbols));
      transform(begin(_symbols), end(_symbols), begin(_arguments), SymEngine::symbol);
      _input.resize(size(_symbols));

      for (auto& setup : _setup)
        setup();

      // set-up final visitor
      auto& visitor = _visitors.emplace_back();
      visitor.init(_arguments,*_se_expression);
      _callbacks.emplace_back([this, &visitor](){
        _result = visitor.call(_input);
      });

    } catch (SymEngine::SymEngineException& e) {
      DUNE_THROW(IOError, e.what());
    }
    _compiled = true;
  }

  [[nodiscard]] RangeField eval() const noexcept override
  {
    assert(_compiled);

    // copy data to input vector
    const std::size_t offset = size(_const_symbols);
    for (std::size_t i = 0; i < size(_input) - offset; ++i)
      _input[i + offset] = *_variables[i];

    // evaluate intermedate expressions
    for (auto& callback : _callbacks)
      callback();

    return _result;
  }

private:
  using Parser::_compiled;
  using Parser::_expression;
  using Parser::_symbols;
  using Parser::_variables;

  std::vector<std::string> _const_symbols;
  mutable std::vector<RangeField> _input;

  SymEngine::RCP<const SymEngine::Basic> _se_expression;
  SymEngine::vec_basic _arguments;

  std::list<SymEngine::LambdaRealDoubleVisitor> _visitors;
  std::vector<std::function<void()>> _setup, _callbacks;
  RangeField _result;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_SYMENGINE_PARSER_HH
