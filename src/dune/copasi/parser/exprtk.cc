#include <dune/copasi/parser/exprtk.hh>

#include <dune/copasi/common/exceptions.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <exprtk.hpp>

#include <array>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#ifndef DUNE_COPASI_EXPRTK_MAX_FUNCTIONS
#error DUNE_COPASI_EXPRTK_MAX_FUNCTIONS is not defined
#endif

namespace Dune::Copasi {

const std::size_t max_functions = DUNE_COPASI_EXPRTK_MAX_FUNCTIONS;

struct ParserData
{
  exprtk::parser<typename ExprTkParser::RangeField> parser;
  exprtk::expression<typename ExprTkParser::RangeField> expression;
  exprtk::symbol_table<typename ExprTkParser::RangeField> symbol_table;
};

// NOLINTBEGIN(altera-struct-pack-align)
template<class F>
struct FunctionID
{
  std::shared_ptr<const void> data;
  std::string name;
  F function;
};
// NOLINTEND(altera-struct-pack-align)

using FunctionID0D = FunctionID<typename ExprTkParser::Function0D>;
using FunctionID1D = FunctionID<typename ExprTkParser::Function1D>;
using FunctionID2D = FunctionID<typename ExprTkParser::Function2D>;
using FunctionID3D = FunctionID<typename ExprTkParser::Function3D>;
using FunctionID4D = FunctionID<typename ExprTkParser::Function4D>;

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
std::array<std::unique_ptr<FunctionID0D>, max_functions> _function0d = {};
std::array<std::unique_ptr<FunctionID1D>, max_functions> _function1d = {};
std::array<std::unique_ptr<FunctionID2D>, max_functions> _function2d = {};
std::array<std::unique_ptr<FunctionID3D>, max_functions> _function3d = {};
std::array<std::unique_ptr<FunctionID4D>, max_functions> _function4d = {};

// recursive mutex are needed because functions may hold parsers that need to be destructed
std::recursive_mutex _mutex = {};
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

ExprTkParser::ExprTkParser()
  : _data{ new ParserData{}, [](void* ptr) { delete static_cast<ParserData*>(ptr); } }
{
}

ExprTkParser::~ExprTkParser()
{
  unregister_functions();
}

void
ExprTkParser::define_constant(const std::string& symbol, const RangeField& value)
{
  static_cast<ParserData*>(_data.get())->symbol_table.add_constant(symbol, value);
};

namespace Impl {

template<std::size_t i>
typename ExprTkParser::RangeField
function_wrapper_0d()
{
  assert(_function0d[i]);
  return _function0d[i]->function();
}

template<std::size_t i>
typename ExprTkParser::RangeField
function_wrapper_1d(typename ExprTkParser::RangeField arg0)
{
  assert(_function1d[i]);
  return _function1d[i]->function(arg0);
}

template<std::size_t i>
typename ExprTkParser::RangeField
function_wrapper_2d(typename ExprTkParser::RangeField arg0, typename ExprTkParser::RangeField arg1)
{
  assert(_function2d[i]);
  return _function2d[i]->function(arg0, arg1);
}

template<std::size_t i>
typename ExprTkParser::RangeField
function_wrapper_3d(typename ExprTkParser::RangeField arg0,
                    typename ExprTkParser::RangeField arg1,
                    typename ExprTkParser::RangeField arg2)
{
  assert(_function3d[i]);
  return _function3d[i]->function(arg0, arg1, arg2);
}

template<std::size_t i>
typename ExprTkParser::RangeField
function_wrapper_4d(typename ExprTkParser::RangeField arg0,
                    typename ExprTkParser::RangeField arg1,
                    typename ExprTkParser::RangeField arg2,
                    typename ExprTkParser::RangeField arg3)
{
  assert(_function4d[i]);
  return _function4d[i]->function(arg0, arg1, arg2, arg3);
}

template<class FunctionID>
void
define_function(auto parser,
                const std::string& symbol,
                auto& function_registry,
                const auto& function)
{
  auto guard = std::unique_lock{ _mutex };
  auto is_assigned = [](const auto& entry) { return bool(entry); };
  // use first slot without assigned functions
  if (auto result = std::ranges::find_if_not(function_registry, is_assigned);
      result != function_registry.end())
    *result = std::make_unique<FunctionID>(FunctionID{ parser, symbol, function });
  else
    throw format_exception(NotImplemented{},
                           "Maximum number of function definitions for ExprTk reached: {}",
                           max_functions);
}

} // namespace Impl

void
ExprTkParser::define_function(const std::string& symbol, const Function0D& function)
{
  Impl::define_function<FunctionID0D>(_data, symbol, _function0d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function1D& function)
{
  Impl::define_function<FunctionID1D>(_data, symbol, _function1d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function2D& function)
{
  Impl::define_function<FunctionID2D>(_data, symbol, _function2d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function3D& function)
{
  Impl::define_function<FunctionID3D>(_data, symbol, _function3d, function);
}

void
ExprTkParser::define_function(const std::string& symbol, const Function4D& function)
{
  Impl::define_function<FunctionID4D>(_data, symbol, _function4d, function);
}

void
ExprTkParser::compile()
{
  assert(not _compiled);
  assert(size(_symbols) == size(_variables));

  ParserData& data = *static_cast<ParserData*>(_data.get());

  // this cast has similar problems as in mu-parser (see comment there)
  // TODO(sospinar): How can we guard this usage of the parser?
  for (std::size_t i = 0; i < size(_symbols); ++i)
    data.symbol_table.add_variable(_symbols[i], *const_cast<double*>(_variables[i]));

  register_functions();

  data.expression.register_symbol_table(data.symbol_table);
  _compiled = data.parser.compile(_expression, data.expression);
  if (!_compiled) {
    throw format_exception(IOError{},
                           "[ExprTk] Failed to compile expression\n"
                           "Expression:    {}\n"
                           "Error message: {}\n",
                           _expression,
                           data.parser.error());
  }
}

[[nodiscard]] ExprTkParser::RangeField
ExprTkParser::operator()() const noexcept
{
  assert(_compiled);
  return static_cast<ParserData*>(_data.get())->expression.value();
}

void
ExprTkParser::register_functions()
{
  ParserData& data = *static_cast<ParserData*>(_data.get());
  auto indices = Dune::range(std::integral_constant<std::size_t, max_functions>{});
  Dune::Hybrid::forEach(indices, [&, guard = std::unique_lock{ _mutex }](auto i) {
    if ((_function0d[i]) and (_function0d[i]->data == _data))
      data.symbol_table.add_function(_function0d[i]->name, Impl::function_wrapper_0d<i>);
    if ((_function1d[i]) and (_function1d[i]->data == _data))
      data.symbol_table.add_function(_function1d[i]->name, Impl::function_wrapper_1d<i>);
    if ((_function2d[i]) and _function2d[i]->data == _data)
      data.symbol_table.add_function(_function2d[i]->name, Impl::function_wrapper_2d<i>);
    if ((_function3d[i]) and _function3d[i]->data == _data)
      data.symbol_table.add_function(_function3d[i]->name, Impl::function_wrapper_3d<i>);
    if ((_function4d[i]) and _function4d[i]->data == _data)
      data.symbol_table.add_function(_function4d[i]->name, Impl::function_wrapper_4d<i>);
  });
}

void
ExprTkParser::unregister_functions()
{
  auto erase_function_from_this = [&](auto& entry) {
    if (entry and entry->data == _data)
      entry.reset();
  };
  auto guard = std::unique_lock{ _mutex };
  std::ranges::for_each(_function0d, erase_function_from_this);
  std::ranges::for_each(_function1d, erase_function_from_this);
  std::ranges::for_each(_function2d, erase_function_from_this);
  std::ranges::for_each(_function3d, erase_function_from_this);
  std::ranges::for_each(_function4d, erase_function_from_this);
}

} // namespace Dune::Copasi
