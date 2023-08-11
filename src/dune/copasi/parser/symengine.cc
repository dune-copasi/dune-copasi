#include <dune/copasi/parser/symengine.hh>
#include <dune/copasi/common/exceptions.hh>

#include <dune/common/exceptions.hh>

#include <symengine/lambda_double.h>
#include <symengine/parser.h>
#include <symengine/parser/sbml/sbml_parser.h>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Dune::Copasi {

SymEngineParser::SymEngineParser(Type parser_type)
  : Parser{}, _parser_type{parser_type}
{
}

void
SymEngineParser::set_expression(const std::string& expression)
{
  Parser::set_expression(expression);
  if (Type::Native == _parser_type)
    _se_expression = SymEngine::parse(_expression);
  else
    _se_expression = SymEngine::parse_sbml(_expression);
}

void
SymEngineParser::define_constant(const std::string& symbol, const RangeField& value)
{
  _input.push_back(value);
  _const_symbols.push_back(symbol);
}

auto
SymEngineParser::setup_function_symbol(const std::string& symbol)
{
  assert(!_se_expression.is_null());
  for (const auto& basic_func : function_symbols(*_se_expression)) {
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
      for (const auto& arg : args) {
        const auto& arg_i = args_i.emplace_back(std::make_shared<std::size_t>(9999999));
        if (has_symbol(*_se_expression, *arg)) {
          _setup.emplace_back([this, arg_i, arg] {
            auto arg_name = SymEngine::rcp_dynamic_cast<const SymEngine::Symbol>(arg)->get_name();
            auto it = find(begin(_symbols), end(_symbols), arg_name);
            if (it != end(_symbols)) {
              (*arg_i) = distance(begin(_symbols), it);
            } else {
              // I am not sure what has_symbol exactly does, so we may get here??
              throw format_exception(
                Exception{}, "Symbol not found in main expression: {}", arg_name);
            }
          });
          continue;
        } else {
          auto arg_name = "__" + fname + "__arg" + std::to_string(count++);
          define_constant(arg_name, 0.);
          *arg_i = size(_input) - 1;
          replace_map[arg] = SymEngine::symbol(arg_name);
          auto& visitor = _visitors.emplace_back();
          _setup.emplace_back([this, &visitor, arg] { visitor.init(_arguments, *arg); });
          _callbacks.emplace_back([&visitor, arg_i, this]() {
            // actual update of argument value
            _input[*arg_i] = visitor.call(_input);
          });
        }
      }
      return args_i;
    }
  }
  throw format_exception(IOError{}, "Function symbol is not contained in expression!");
}

void
SymEngineParser::define_function(const std::string& symbol, const Function0D& function)
{
  auto args_i = setup_function_symbol(symbol);

  if (size(args_i) != 1)
    throw format_exception(IOError{}, "Function arguments do not match with defined function");

  _callbacks.emplace_back([args_i, this, function = function]() {
    // update function return value
    _input[*args_i[0]] = function();
  });
}

void
SymEngineParser::define_function(const std::string& symbol, const Function1D& function)
{
  auto args_i = setup_function_symbol(symbol);

  if (size(args_i) != 2)
    throw format_exception(IOError{}, "Function arguments do not match with defined function");

  _callbacks.emplace_back([args_i, this, function = function]() {
    // update function return value
    _input[*args_i[0]] = function(_input[*args_i[1]]);
  });
}

void
SymEngineParser::define_function(const std::string& symbol, const Function2D& function)
{
  auto args_i = setup_function_symbol(symbol);

  if (size(args_i) != 3)
    throw format_exception(IOError{}, "Function arguments do not match with defined function");

  _callbacks.emplace_back([args_i, this, function = function]() {
    // update function return value
    _input[*args_i[0]] = function(_input[*args_i[1]], _input[*args_i[2]]);
  });
}

void
SymEngineParser::define_function(const std::string& symbol, const Function3D& function)
{
  auto args_i = setup_function_symbol(symbol);

  if (size(args_i) != 4)
    throw format_exception(IOError{}, "Function arguments do not match with defined function");

  _callbacks.emplace_back([args_i, this, function = function]() {
    // update function return value
    _input[*args_i[0]] = function(_input[*args_i[1]], _input[*args_i[2]], _input[*args_i[3]]);
  });
}

void
SymEngineParser::define_function(const std::string& symbol, const Function4D& function)
{
  auto args_i = setup_function_symbol(symbol);

  if (size(args_i) != 5)
    throw format_exception(IOError{}, "Function arguments do not match with defined function");

  _callbacks.emplace_back([args_i, this, function = function]() {
    // update function return value
    _input[*args_i[0]] = function(_input[*args_i[1]], _input[*args_i[2]], _input[*args_i[3]], _input[*args_i[4]]);
  });
}

void
SymEngineParser::compile()
{
  Parser::compile();
  assert(not _compiled);
  assert(size(_symbols) == size(_variables));

  try {
    _symbols.insert(begin(_symbols), begin(_const_symbols), end(_const_symbols));
    _arguments.resize(size(_symbols));
    transform(begin(_symbols), end(_symbols), begin(_arguments), SymEngine::symbol);
    _input.resize(size(_symbols));

    for (const auto& setup : _setup)
      setup();

    // set-up final visitor
    auto& visitor = _visitors.emplace_back();
    visitor.init(_arguments, *_se_expression);
    _callbacks.emplace_back([this, &visitor]() { _result = visitor.call(_input); });

  } catch (SymEngine::SymEngineException& e) {
    throw format_exception(IOError{}, "{}", e.what());
  }
  _compiled = true;
}

[[nodiscard]] SymEngineParser::RangeField
SymEngineParser::operator()() const noexcept
{
  assert(_compiled);

  // copy data to input vector
  const std::size_t offset = size(_const_symbols);
  for (std::size_t i = 0; i < size(_input) - offset; ++i)
    _input[i + offset] = *_variables[i];

  // evaluate intermedate expressions
  for (const auto& callback : _callbacks)
    callback();

  return _result;
}

} // namespace Dune::Copasi
