#ifndef DUNE_COPASI_PARSER_PARSER_HH
#define DUNE_COPASI_PARSER_PARSER_HH

#include <dune-copasi-config.h>

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

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

  using Function0D = std::function<RangeField()>;
  using Function1D = std::function<RangeField(RangeField)>;
  using Function2D = std::function<RangeField(RangeField, RangeField)>;
  using Function3D = std::function<RangeField(RangeField, RangeField, RangeField)>;
  using Function4D = std::function<RangeField(RangeField, RangeField, RangeField, RangeField)>;

  virtual void set_expression(const std::string& expression);

  [[nodiscard]] std::string expression() const;

  virtual void define_variable(const std::string& symbol, RangeField const* value);

  virtual void define_constant(const std::string& symbol, const RangeField& value) = 0;

  virtual void define_function(const std::string& symbol, const Function0D& function) = 0;
  virtual void define_function(const std::string& symbol, const Function1D& function) = 0;
  virtual void define_function(const std::string& symbol, const Function2D& function) = 0;
  virtual void define_function(const std::string& symbol, const Function3D& function) = 0;
  virtual void define_function(const std::string& symbol, const Function4D& function) = 0;

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
