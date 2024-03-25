#ifndef DUNE_COPASI_MU_PARSER_HH
#define DUNE_COPASI_MU_PARSER_HH

#include <dune/copasi/parser/parser.hh>

#include <memory>
#include <string>

namespace Dune::Copasi {

class MuParser final : public Parser
{
public:
  MuParser();

  MuParser(const MuParser&) = delete;
  MuParser(MuParser&&) = default;

  MuParser& operator=(const MuParser&) = delete;
  MuParser& operator=(MuParser&&) = default;

  ~MuParser() override final;

  using RangeField = typename Parser::RangeField;

  using Function0D = typename Parser::Function0D;
  using Function1D = typename Parser::Function1D;
  using Function2D = typename Parser::Function2D;
  using Function3D = typename Parser::Function3D;

  void define_constant(const std::string& symbol, const RangeField& value) final;
  void define_function(const std::string& symbol, Function0D&& function) final;
  void define_function(const std::string& symbol, Function1D&& function) final;
  void define_function(const std::string& symbol, Function2D&& function) final;
  void define_function(const std::string& symbol, Function3D&& function) final;
  void define_function(const std::string& symbol, Function4D&& function) final;

  void compile() final;

  [[nodiscard]] RangeField operator()() const noexcept final;

private:
  void register_functions();

  void unregister_functions();
  using Parser::_compiled;
  using Parser::_expression;
  using Parser::_symbols;
  using Parser::_variables;

  std::shared_ptr<void> _parser;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MU_PARSER_HH
