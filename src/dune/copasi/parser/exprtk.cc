#include <dune/copasi/parser/exprtk.hh>

#include <dune/copasi/common/exceptions.hh>

#include <exprtk.hpp>

#include <map>
#include <mutex>
#include <string>
#include <vector>

#ifndef DUNE_COPASI_EXPRTK_MAX_FUNCTIONS
#define DUNE_COPASI_EXPRTK_MAX_FUNCTIONS 500
#endif


namespace Dune::Copasi {

const std::size_t max_functions = DUNE_COPASI_EXPRTK_MAX_FUNCTIONS;

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

namespace Impl {

template<class FunctionID>
void define_function(auto parser, const std::string& symbol, auto& function_registry, const auto& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i) {
    if (not function_registry[i]) {
      function_registry[i] = std::make_unique<FunctionID>(FunctionID{ parser, symbol, function });
      return;
    }
  }
  throw format_exception(
    NotImplemented{}, "Maximum number of functions reached: {}", max_functions);
}

} // namespace Impl


void
ExprTkParser::define_function(const std::string& symbol, const Function0D& function)
{
  Impl::define_function<FunctionID0D>(_parser, symbol, _function0d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function1D& function)
{
  Impl::define_function<FunctionID1D>(_parser, symbol, _function1d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function2D& function)
{
  Impl::define_function<FunctionID2D>(_parser, symbol, _function2d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function3D& function)
{
  Impl::define_function<FunctionID3D>(_parser, symbol, _function3d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function4D& function)
{
  Impl::define_function<FunctionID4D>(_parser, symbol, _function4d, function);
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
    auto i4 = _function4d.find(i);

    if (i0 != end(_function0d) and (i0->second) and (i0->second->parser == _parser))
      _symbol_table.add_function(i0->second->name, function_wrapper_0d<i>);
    if (i1 != end(_function1d) and (i1->second) and (i1->second->parser == _parser))
      _symbol_table.add_function(i1->second->name, function_wrapper_1d<i>);
    if (i2 != end(_function2d) and (i2->second) and i2->second->parser == _parser)
      _symbol_table.add_function(i2->second->name, function_wrapper_2d<i>);
    if (i3 != end(_function3d) and (i3->second) and i3->second->parser == _parser)
      _symbol_table.add_function(i3->second->name, function_wrapper_3d<i>);
    if (i4 != end(_function4d) and (i4->second) and i4->second->parser == _parser)
      _symbol_table.add_function(i4->second->name, function_wrapper_4d<i>);
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
    auto i4 = _function4d.find(i);

    if (i0 != end(_function0d) and (i0->second) and (i0->second->parser == _parser))
      _function0d.erase(i0);
    if (i1 != end(_function1d) and (i1->second) and (i1->second->parser == _parser))
      _function1d.erase(i1);
    if (i2 != end(_function2d) and (i2->second) and i2->second->parser == _parser)
      _function2d.erase(i2);
    if (i3 != end(_function3d) and (i3->second) and i3->second->parser == _parser)
      _function3d.erase(i3);
    if (i4 != end(_function4d) and (i4->second) and i4->second->parser == _parser)
      _function4d.erase(i4);
  });
}

} // namespace Dune::Copasi
