#ifndef DUNE_COPASI_EXPRTK_PARSER_HH
#define DUNE_COPASI_EXPRTK_PARSER_HH

#include <dune/copasi/parser/parser.hh>

#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <memory>

namespace Dune::Copasi {

class ExprTkParser final : public Parser
{
  static const std::size_t max_functions = 100;

public:
  ExprTkParser();

  ExprTkParser(const ExprTkParser&) = delete;
  ExprTkParser(ExprTkParser&&) = default;

  ExprTkParser& operator=(const ExprTkParser&) = delete;
  ExprTkParser& operator=(ExprTkParser&&) = default;

  ~ExprTkParser() override final;

  using RangeField = typename Parser::RangeField;

  using Function0D = typename Parser::Function0D;
  using Function1D = typename Parser::Function1D;
  using Function2D = typename Parser::Function2D;
  using Function3D = typename Parser::Function3D;
  using Function4D = typename Parser::Function4D;

  void define_constant(const std::string& symbol, const RangeField& value) override final;

  void define_function(const std::string& symbol, const Function0D& function) override final;

  void define_function(const std::string& symbol, const Function1D& function) override final;

  void define_function(const std::string& symbol, const Function2D& function) override final;

  void define_function(const std::string& symbol, const Function3D& function) override final;

  void define_function(const std::string& symbol, const Function4D& function) override final;

  void compile() override final;

  [[nodiscard]] RangeField operator()() const noexcept override final;

private:
  void register_functions();
  void unregister_functions();

  using Parser::_compiled;
  using Parser::_expression;
  using Parser::_symbols;
  using Parser::_variables;

  std::shared_ptr<void> _data;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_EXPRTK_PARSER_HH
