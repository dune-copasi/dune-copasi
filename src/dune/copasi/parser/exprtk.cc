#include <dune/copasi/parser/exprtk.hh>

#include <exprtk.hpp>

#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace Dune::Copasi {

ExprTkParser::ExprTkParser()
  : Parser{}
  , _parser{ std::make_shared<exprtk::parser<RangeField>>() }
{
}

ExprTkParser::~ExprTkParser()
{
  unregister_functions();
}

void
ExprTkParser::define_constant(const std::string& symbol, const RangeField& value)
{
  _symbol_table.add_constant(symbol, value);
};

void
ExprTkParser::define_function(const std::string& symbol, const Function0D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i)
    if (not _function0d[i]) {
      _function0d[i] = std::make_unique<FunctionID0D>(FunctionID0D{ _parser, symbol, function });
      return;
    }
  DUNE_THROW(Exception, "\tNot enough instantiated functions");
}

void
ExprTkParser::define_function(const std::string& symbol, const Function1D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i)
    if (not _function1d[i]) {
      _function1d[i] = std::make_unique<FunctionID1D>(FunctionID1D{ _parser, symbol, function });
      return;
    }
  DUNE_THROW(Exception, "\tNot enough instantiated functions");
}

void
ExprTkParser::define_function(const std::string& symbol, const Function2D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i)
    if (not _function2d[i]) {
      _function2d[i] = std::make_unique<FunctionID2D>(FunctionID2D{ _parser, symbol, function });
      return;
    }
  DUNE_THROW(Exception, "\tNot enough instantiated functions");
}

void
ExprTkParser::define_function(const std::string& symbol, const Function3D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i)
    if (not _function3d[i]) {
      _function3d[i] = std::make_unique<FunctionID3D>(FunctionID3D{ _parser, symbol, function });
      return;
    }
  DUNE_THROW(Exception, "\tNot enough instantiated functions");
}

void
ExprTkParser::compile()
{
  Parser::compile();
  assert(not _compiled);
  assert(size(_symbols) == size(_variables));

  for (std::size_t i = 0; i < size(_symbols); ++i)
    _symbol_table.add_variable(_symbols[i], *_variables[i]);

  register_functions();

  _tk_expression.register_symbol_table(_symbol_table);
  _parser->compile(_expression, _tk_expression);
  _compiled = true;
}

[[nodiscard]] MuParser::RangeField
ExprTkParser::operator()() const noexcept
{
  assert(_compiled);
  return _tk_expression.value();
}

void
ExprTkParser::register_functions()
{
  auto indices = Dune::range(std::integral_constant<std::size_t, max_functions>{});
  Dune::Hybrid::forEach(indices, [&, guard = std::unique_lock{ _mutex }](auto i) {
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
}

void
ExprTkParser::unregister_functions()
{
  auto indices = Dune::range(std::integral_constant<std::size_t, max_functions>{});
  Dune::Hybrid::forEach(indices, [&, guard = std::unique_lock{ _mutex }](auto i) {
    auto i0 = _function0d.find(i);
    auto i1 = _function1d.find(i);
    auto i2 = _function2d.find(i);
    auto i3 = _function3d.find(i);

    if (i0 != end(_function0d) and (i0->second) and (i0->second->parser == _parser))
      _function0d.erase(i0);
    if (i1 != end(_function1d) and (i1->second) and (i1->second->parser == _parser))
      _function1d.erase(i1);
    if (i2 != end(_function2d) and (i2->second) and i2->second->parser == _parser)
      _function2d.erase(i2);
    if (i3 != end(_function3d) and (i3->second) and i3->second->parser == _parser)
      _function3d.erase(i3);
  });
}

} // namespace Dune::Copasi
