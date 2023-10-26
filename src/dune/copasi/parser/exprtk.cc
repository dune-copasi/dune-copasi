#include <dune/copasi/parser/exprtk.hh>

#include <dune/copasi/common/exceptions.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <exprtk.hpp>

#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <memory>

#ifndef DUNE_COPASI_EXPRTK_MAX_FUNCTIONS
#define DUNE_COPASI_EXPRTK_MAX_FUNCTIONS 500
#endif


namespace Dune::Copasi {

const std::size_t max_functions = DUNE_COPASI_EXPRTK_MAX_FUNCTIONS;

struct ParserData {
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
std::map<std::size_t, std::unique_ptr<FunctionID0D>> _function0d = {};
std::map<std::size_t, std::unique_ptr<FunctionID1D>> _function1d = {};
std::map<std::size_t, std::unique_ptr<FunctionID2D>> _function2d = {};
std::map<std::size_t, std::unique_ptr<FunctionID3D>> _function3d = {};
std::map<std::size_t, std::unique_ptr<FunctionID4D>> _function4d = {};

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
  data.parser.compile(_expression, data.expression);
  _compiled = true;
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
    auto i0 = _function0d.find(i);
    auto i1 = _function1d.find(i);
    auto i2 = _function2d.find(i);
    auto i3 = _function3d.find(i);
    auto i4 = _function4d.find(i);

    if (i0 != end(_function0d) and (i0->second) and (i0->second->data == _data))
      data.symbol_table.add_function(i0->second->name, Impl::function_wrapper_0d<i>);
    if (i1 != end(_function1d) and (i1->second) and (i1->second->data == _data))
      data.symbol_table.add_function(i1->second->name, Impl::function_wrapper_1d<i>);
    if (i2 != end(_function2d) and (i2->second) and i2->second->data == _data)
      data.symbol_table.add_function(i2->second->name, Impl::function_wrapper_2d<i>);
    if (i3 != end(_function3d) and (i3->second) and i3->second->data == _data)
      data.symbol_table.add_function(i3->second->name, Impl::function_wrapper_3d<i>);
    if (i4 != end(_function4d) and (i4->second) and i4->second->data == _data)
      data.symbol_table.add_function(i4->second->name, Impl::function_wrapper_4d<i>);
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

    if (i0 != end(_function0d) and (i0->second) and (i0->second->data == _data))
      _function0d.erase(i0);
    if (i1 != end(_function1d) and (i1->second) and (i1->second->data == _data))
      _function1d.erase(i1);
    if (i2 != end(_function2d) and (i2->second) and i2->second->data == _data)
      _function2d.erase(i2);
    if (i3 != end(_function3d) and (i3->second) and i3->second->data == _data)
      _function3d.erase(i3);
    if (i4 != end(_function4d) and (i4->second) and i4->second->data == _data)
      _function4d.erase(i4);
  });
}

} // namespace Dune::Copasi
