#include <algorithm>
#include <dune/copasi/parser/mu.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/parser/parser.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <spdlog/spdlog.h>

#include <muParser.h>

#include <array>
#include <cassert>
#include <cstddef>
#include <exception>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

#ifndef DUNE_COPASI_MUPARSER_MAX_FUNCTIONS
#error DUNE_COPASI_MUPARSER_MAX_FUNCTIONS is not defined
#endif

namespace Dune::Copasi {

namespace {

const std::size_t max_functions = DUNE_COPASI_MUPARSER_MAX_FUNCTIONS;

// NOLINTBEGIN(altera-struct-pack-align)
template<class F>
struct FunctionID
{
  std::shared_ptr<const void> parser = nullptr;
  std::string name = {};
  F function = {};
};
// NOLINTEND(altera-struct-pack-align)

using FunctionID0D = FunctionID<typename MuParser::Function0D>;
using FunctionID1D = FunctionID<typename MuParser::Function1D>;
using FunctionID2D = FunctionID<typename MuParser::Function2D>;
using FunctionID3D = FunctionID<typename MuParser::Function3D>;
using FunctionID4D = FunctionID<typename MuParser::Function4D>;

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
std::array<FunctionID0D, max_functions> _function0d = {};
std::array<FunctionID1D, max_functions> _function1d = {};
std::array<FunctionID2D, max_functions> _function2d = {};
std::array<FunctionID3D, max_functions> _function3d = {};
std::array<FunctionID4D, max_functions> _function4d = {};

// recursive mutex are needed because functions may hold parsers that need to be destructed
std::recursive_mutex _mutex = {};
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_0d()
{
  return _function0d[i].function();
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_1d(typename MuParser::RangeField arg0)
{
  return _function1d[i].function(arg0);
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_2d(typename MuParser::RangeField arg0, typename MuParser::RangeField arg1)
{
  return _function2d[i].function(arg0, arg1);
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_3d(typename MuParser::RangeField arg0,
                    typename MuParser::RangeField arg1,
                    typename MuParser::RangeField arg2)
{
  return _function3d[i].function(arg0, arg1, arg2);
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_4d(typename MuParser::RangeField arg0,
                    typename MuParser::RangeField arg1,
                    typename MuParser::RangeField arg2,
                    typename MuParser::RangeField arg3)
{
  return _function4d[i].function(arg0, arg1, arg2, arg3);
}

} // namespace

MuParser::MuParser()
  : _parser{ new mu::Parser{}, [](void* ptr) { delete static_cast<mu::Parser*>(ptr); } }
{
}

MuParser::~MuParser()
{
  unregister_functions();
}

void
MuParser::define_constant(const std::string& symbol, const RangeField& value)
{
  static_cast<mu::Parser*>(_parser.get())->DefineConst(symbol, value);
}

namespace Impl {

template<class F>
void
define_function(std::shared_ptr<const void> parser,
                const std::string& symbol,
                auto& function_registry,
                F&& function)
{
  auto guard = std::unique_lock{ _mutex };
  auto is_assigned = [](const auto& entry) { return bool(entry.parser); };
  // use first slot without assigned functions
  if (auto result = std::ranges::find_if_not(function_registry, is_assigned);
      result != function_registry.end())
    *result = FunctionID<F>{ std::move(parser), symbol, std::forward<F>(function) };
  else
    throw format_exception(NotImplemented{},
                           "Maximum number of function definitions for muParser reached: {}",
                           max_functions);
}

} // namespace Impl

void
MuParser::define_function(const std::string& symbol, Function0D&& function)
{
  Impl::define_function(_parser, symbol, _function0d, std::move(function));
}

void
MuParser::define_function(const std::string& symbol, Function1D&& function)
{
  Impl::define_function(_parser, symbol, _function1d, std::move(function));
}

void
MuParser::define_function(const std::string& symbol, Function2D&& function)
{
  Impl::define_function(_parser, symbol, _function2d, std::move(function));
}

void
MuParser::define_function(const std::string& symbol, Function3D&& function)
{
  Impl::define_function(_parser, symbol, _function3d, std::move(function));
}

void
MuParser::define_function(const std::string& symbol, Function4D&& function)
{
  Impl::define_function(_parser, symbol, _function4d, std::move(function));
}

void
MuParser::compile()
{
  assert(not _compiled);
  assert(size(_symbols) == size(_variables));

  auto& mu_parser = *static_cast<mu::Parser*>(_parser.get());

  // notice that we are casting away the constness of the variable because
  // that's muParser signature. In practice, this means that the muParser
  // expression should not attempt to modify the variable at all (e.g., var=10),
  // otherwise, it will result in undefined behavior!
  // TODO(sospinar): How can we guard this usage of the parser?
  for (std::size_t i = 0; i < size(_symbols); ++i) {
    mu_parser.DefineVar(
      _symbols[i],
      const_cast<double*>(_variables[i])); // NOLINT(cppcoreguidelines-pro-type-const-cast)
  }

  register_functions();

  mu_parser.EnableOptimizer();
  mu_parser.SetExpr(_expression);
  _compiled = true;

  try {
    mu_parser.Eval();
  } catch (mu::Parser::exception_type& e) {
    throw format_exception(IOError{},
                           "Evaluating muParser expression failed:\n"
                           "  Parsed expression:   {}\n"
                           "  Token:               {}\n"
                           "  Error position:      {}\n"
                           "  Error code:          {}\n"
                           "  Error message:       {}\n",
                           e.GetExpr(),
                           e.GetToken(),
                           e.GetPos(),
                           static_cast<int>(e.GetCode()),
                           e.GetMsg());
  }
}

[[nodiscard]] MuParser::RangeField
MuParser::operator()() const noexcept
{
  try {
    assert(_compiled);
    return static_cast<mu::Parser*>(_parser.get())->Eval();
  } catch (mu::Parser::exception_type& e) {
    spdlog::error("Evaluating muParser expression failed:\n"
                  "  Parsed expression:   {}\n"
                  "  Token:               {}\n"
                  "  Error position:      {}\n"
                  "  Error code:          {}\n"
                  "  Error message:       {}\n",
                  e.GetExpr(),
                  e.GetToken(),
                  e.GetPos(),
                  static_cast<int>(e.GetCode()),
                  e.GetMsg());
    std::terminate();
  } catch (...) {
    std::terminate();
  }
}

void
MuParser::register_functions()
{
  auto indices = std::make_index_sequence<max_functions>{};
  auto& mu_parser = *static_cast<mu::Parser*>(_parser.get());
  auto guard = std::unique_lock{ _mutex };
  Dune::Hybrid::forEach(indices, [&](auto idx) {
    if (_function0d[idx].parser == _parser)
      mu_parser.DefineFun(_function0d[idx].name, function_wrapper_0d<idx>);
    if (_function1d[idx].parser == _parser)
      mu_parser.DefineFun(_function1d[idx].name, function_wrapper_1d<idx>);
    if (_function2d[idx].parser == _parser)
      mu_parser.DefineFun(_function2d[idx].name, function_wrapper_2d<idx>);
    if (_function3d[idx].parser == _parser)
      mu_parser.DefineFun(_function3d[idx].name, function_wrapper_3d<idx>);
    if (_function4d[idx].parser == _parser)
      mu_parser.DefineFun(_function4d[idx].name, function_wrapper_4d<idx>);
  });
}

void
MuParser::unregister_functions()
{
  auto erase_function_from_this = [&](auto& entry) {
    if (entry.parser == _parser)
      entry = {};
  };
  auto guard = std::unique_lock{ _mutex };
  std::ranges::for_each(_function0d, erase_function_from_this);
  std::ranges::for_each(_function1d, erase_function_from_this);
  std::ranges::for_each(_function2d, erase_function_from_this);
  std::ranges::for_each(_function3d, erase_function_from_this);
  std::ranges::for_each(_function4d, erase_function_from_this);
}

} // namespace Dune::Copasi
