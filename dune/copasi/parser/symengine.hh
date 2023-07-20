#ifndef DUNE_COPASI_SYMENGINE_PARSER_HH
#define DUNE_COPASI_SYMENGINE_PARSER_HH

#include <dune/copasi/parser/parser.hh>

#include <symengine/lambda_double.h>
#include <symengine/parser.h>
#include <symengine/parser/sbml/sbml_parser.h>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Dune::Copasi {

class SymEngineParser final : public Parser
{

public:
  enum class Type
  {
    Native,
    SBML
  };

  using RangeField = typename Parser::RangeField;

  using Function0D = typename Parser::Function0D;
  using Function1D = typename Parser::Function1D;
  using Function2D = typename Parser::Function2D;
  using Function3D = typename Parser::Function3D;

  SymEngineParser(Type parser_type = Type::Native);

  void set_expression(const std::string& expression) override final;

  void define_constant(const std::string& symbol, const RangeField& value) override final;

  auto setup_function_symbol(const std::string& symbol);

  void define_function(const std::string& symbol, const Function0D& function) override final;
  void define_function(const std::string& symbol, const Function1D& function) override final;
  void define_function(const std::string& symbol, const Function2D& function) override final;
  void define_function(const std::string& symbol, const Function3D& function) override final;

  void compile() override final;

  [[nodiscard]] RangeField operator()() const noexcept override final;

private:
  using Parser::_compiled;
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