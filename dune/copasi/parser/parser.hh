#ifndef DUNE_COPASI_PARSER_PARSER_HH
#define DUNE_COPASI_PARSER_PARSER_HH

#include <dune-copasi-config.hh>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include <function2/function2.hpp>

namespace Dune::Copasi {

class Parser
{

public:
  Parser() = default;

  Parser(const Parser&) = delete;
  Parser(Parser&&) = default;

  Parser& operator=(const Parser&) = delete;
  Parser& operator=(Parser&&) = default;

  virtual ~Parser() = default;

  using RangeField = double;

  using Function0D = fu2::unique_function<RangeField() const>;
  using Function1D = fu2::unique_function<RangeField(RangeField) const>;
  using Function2D = fu2::unique_function<RangeField(RangeField, RangeField) const>;
  using Function3D = fu2::unique_function<RangeField(RangeField, RangeField, RangeField) const>;
  using Function4D = fu2::unique_function<RangeField(RangeField, RangeField, RangeField, RangeField) const>;

  virtual void set_expression(const std::string& expression);

  [[nodiscard]] std::string expression() const;

  virtual void define_variable(const std::string& symbol, RangeField const* value);

  virtual void define_constant(const std::string& symbol, const RangeField& value) = 0;

  virtual void define_function(const std::string& symbol, Function0D&& function) = 0;
  virtual void define_function(const std::string& symbol, Function1D&& function) = 0;
  virtual void define_function(const std::string& symbol, Function2D&& function) = 0;
  virtual void define_function(const std::string& symbol, Function3D&& function) = 0;
  virtual void define_function(const std::string& symbol, Function4D&& function) = 0;

  [[nodiscard]] bool compiled() const;

  virtual void compile() = 0;

  [[nodiscard]] virtual RangeField operator()() const noexcept = 0;

protected:
  bool _compiled = false;
  std::string _expression;
  std::vector<std::string> _symbols;
  std::vector<RangeField const*> _variables;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_PARSER_HH
