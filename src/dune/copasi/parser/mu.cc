#include <dune/copasi/parser/mu.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/parser/parser.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

#include <spdlog/spdlog.h>

#include <muParser.h>

#include <cassert>
#include <cstddef>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

namespace Dune::Copasi {

namespace {

const std::size_t max_functions = 100;

// NOLINTBEGIN(altera-struct-pack-align)
template<class F>
struct FunctionID
{
  std::shared_ptr<const void> parser;
  std::string name;
  F function;
};
// NOLINTEND(altera-struct-pack-align)

using FunctionID0D = FunctionID<typename MuParser::Function0D>;
using FunctionID1D = FunctionID<typename MuParser::Function1D>;
using FunctionID2D = FunctionID<typename MuParser::Function2D>;
using FunctionID3D = FunctionID<typename MuParser::Function3D>;

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
std::map<std::size_t, std::unique_ptr<FunctionID0D>> _function0d = {};
std::map<std::size_t, std::unique_ptr<FunctionID1D>> _function1d = {};
std::map<std::size_t, std::unique_ptr<FunctionID2D>> _function2d = {};
std::map<std::size_t, std::unique_ptr<FunctionID3D>> _function3d = {};

// recursive mutex are needed because functions may hold parsers that need to be destructed
std::recursive_mutex _mutex = {};
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_0d()
{
  assert(_function0d[i]);
  return _function0d[i]->function();
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_1d(typename MuParser::RangeField arg0)
{
  assert(_function1d[i]);
  return _function1d[i]->function(arg0);
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_2d(typename MuParser::RangeField arg0, typename MuParser::RangeField arg1)
{
  assert(_function2d[i]);
  return _function2d[i]->function(arg0, arg1);
}

template<std::size_t i>
typename MuParser::RangeField
function_wrapper_3d(typename MuParser::RangeField arg0,
                    typename MuParser::RangeField arg1,
                    typename MuParser::RangeField arg2)
{
  assert(_function3d[i]);
  return _function3d[i]->function(arg0, arg1, arg2);
}

} // namespace

MuParser::MuParser()
  : _parser{ new mu::Parser{}, [](void* ptr) { delete reinterpret_cast<mu::Parser*>(ptr); } }
{
}

MuParser::~MuParser()
{
  unregister_functions();
}

void
MuParser::define_constant(const std::string& symbol, const RangeField& value)
{
  reinterpret_cast<mu::Parser*>(_parser.get())->DefineConst(symbol, value);
}

void
MuParser::define_function(const std::string& symbol, const Function0D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i) {
    if (not _function0d[i]) {
      _function0d[i] = std::make_unique<FunctionID0D>(FunctionID0D{ _parser, symbol, function });
      return;
    }
  }
  throw format_exception(
    NotImplemented{}, "Maximum number of functions reached: {}", max_functions);
}

void
MuParser::define_function(const std::string& symbol, const Function1D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i) {
    if (not _function1d[i]) {
      _function1d[i] = std::make_unique<FunctionID1D>(FunctionID1D{ _parser, symbol, function });
      return;
    }
  }
  throw format_exception(
    NotImplemented{}, "Maximum number of functions reached: {}", max_functions);
}

void
MuParser::define_function(const std::string& symbol, const Function2D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i) {
    if (not _function2d[i]) {
      _function2d[i] = std::make_unique<FunctionID2D>(FunctionID2D{ _parser, symbol, function });
      return;
    }
  }
  throw format_exception(
    NotImplemented{}, "Maximum number of functions reached: {}", max_functions);
}

void
MuParser::define_function(const std::string& symbol, const Function3D& function)
{
  auto guard = std::unique_lock{ _mutex };
  for (std::size_t i = 0; i < max_functions; ++i) {
    if (not _function3d[i]) {
      _function3d[i] = std::make_unique<FunctionID3D>(FunctionID3D{ _parser, symbol, function });
      return;
    }
  }
  throw format_exception(
    NotImplemented{}, "Maximum number of functions reached: {}", max_functions);
}

void
MuParser::compile()
{
  assert(not _compiled);
  assert(size(_symbols) == size(_variables));

  auto& mu_parser = *reinterpret_cast<mu::Parser*>(_parser.get());

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
    return reinterpret_cast<mu::Parser*>(_parser.get())->Eval();
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
  auto& mu_parser = *reinterpret_cast<mu::Parser*>(_parser.get());
  Dune::Hybrid::forEach(indices, [&, guard = std::unique_lock{ _mutex }](auto idx) {
    auto idx_0d = _function0d.find(idx);
    auto idx_1d = _function1d.find(idx);
    auto idx_2d = _function2d.find(idx);
    auto idx_3d = _function3d.find(idx);

    if (idx_0d != end(_function0d) and (idx_0d->second) and (idx_0d->second->parser == _parser)) {
      mu_parser.DefineFun(idx_0d->second->name, function_wrapper_0d<idx>);
    }
    if (idx_1d != end(_function1d) and (idx_1d->second) and (idx_1d->second->parser == _parser)) {
      mu_parser.DefineFun(idx_1d->second->name, function_wrapper_1d<idx>);
    }
    if (idx_2d != end(_function2d) and (idx_2d->second) and idx_2d->second->parser == _parser) {
      mu_parser.DefineFun(idx_2d->second->name, function_wrapper_2d<idx>);
    }
    if (idx_3d != end(_function3d) and (idx_3d->second) and idx_3d->second->parser == _parser) {
      mu_parser.DefineFun(idx_3d->second->name, function_wrapper_3d<idx>);
    }
  });
}

void
MuParser::unregister_functions()
{
  auto indices = std::make_index_sequence<max_functions>{};
  Dune::Hybrid::forEach(indices, [&, guard = std::unique_lock{ _mutex }](auto idx) {
    auto idx_0d = _function0d.find(idx);
    auto idx_1d = _function1d.find(idx);
    auto idx_2d = _function2d.find(idx);
    auto idx_3d = _function3d.find(idx);

    if (idx_0d != end(_function0d) and (idx_0d->second) and (idx_0d->second->parser == _parser)) {
      _function0d.erase(idx_0d);
    }
    if (idx_1d != end(_function1d) and (idx_1d->second) and (idx_1d->second->parser == _parser)) {
      _function1d.erase(idx_1d);
    }
    if (idx_2d != end(_function2d) and (idx_2d->second) and idx_2d->second->parser == _parser) {
      _function2d.erase(idx_2d);
    }
    if (idx_3d != end(_function3d) and (idx_3d->second) and idx_3d->second->parser == _parser) {
      _function3d.erase(idx_3d);
    }
  });
}

} // namespace Dune::Copasi
