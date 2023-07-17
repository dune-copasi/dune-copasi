#include <dune/copasi/parser/parser.hh>

#include <dune/common/exceptions.hh>

#include <fmt/core.h>
#include <fmt/format.h>

#include <algorithm>
#include <string>

namespace Dune::Copasi {

void
Parser::set_expression(const std::string& expression)
{
  _expression = expression;
}

[[nodiscard]] auto
Parser::expression() const -> std::string
{
  return _expression;
}

void
Parser::define_variable(const std::string& symbol, RangeField const* value)
{
  auto symb_it = std::find(_symbols.begin(), _symbols.end(), symbol);
  if (symb_it != _symbols.end()) {
    DUNE_THROW(Exception, fmt::format("\tSymbol '{}' is already defined", symbol));
  }
  _symbols.push_back(symbol);
  _variables.emplace_back(value);
}

} // namespace Dune::Copasi
