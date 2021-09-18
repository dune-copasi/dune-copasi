#ifndef DUNE_COPASI_PARSER_HH
#define DUNE_COPASI_PARSER_HH

// #include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>

#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Dune::Copasi {

class Parser
{

public:
  Parser(const std::string& expression)
    : _compiled{ false }
    , _expression{ expression }
  {}

  virtual ~Parser(){};

  using RangeField = double;

  using Function0D = std::function<RangeField()>;
  using Function1D = std::function<RangeField(RangeField)>;
  using Function2D = std::function<RangeField(RangeField, RangeField)>;
  using Function3D =
    std::function<RangeField(RangeField, RangeField, RangeField)>;

  virtual void define_variable(const std::string& symbol, RangeField* value)
  {
    auto it = find(begin(_symbols), end(_symbols), symbol);
    if (it != end(_symbols))
      DUNE_THROW(Exception, "");
    _symbols.push_back(symbol);
    _variables.emplace_back(value);
  }

  virtual void define_constant(const std::string& symbol,
                               const RangeField& value) = 0;

  virtual void define_function(const std::string& symbol,
                               Function0D&& function) = 0;
  virtual void define_function(const std::string& symbol,
                               Function1D&& function) = 0;
  virtual void define_function(const std::string& symbol,
                               Function2D&& function) = 0;
  virtual void define_function(const std::string& symbol,
                               Function3D&& function) = 0;

  [[nodiscard]] bool compiled() const { return _compiled; }

  virtual void compile() {}

  [[nodiscard]] virtual RangeField eval() const noexcept = 0;

protected:
  bool _compiled;
  std::string _expression;
  std::vector<std::string> _symbols;
  std::vector<RangeField*> _variables;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_HH
