#include <dune-copasi-config.hh>

#include <dune/copasi/parser/exprtk.hh>

#include <dune/copasi/common/exceptions.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>

#include <exprtk.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#if (__cpp_lib_memory_resource >= 201603L) && (__cpp_lib_polymorphic_allocator >= 201902L)
#define DUNE_COPASI_PARSER_EXPRTK_USES_MEMORY_RESOURCE
#include <memory_resource>
#endif

namespace Dune::Copasi {

struct ParserData
{
#ifdef DUNE_COPASI_PARSER_EXPRTK_USES_MEMORY_RESOURCE
  // Important: this memory resource will hold the contents of the smart pointers of 'functions' contiguously.
  // Hence, this memory resource needs to be destroyed **after** 'functions' is completely destroyed.
  // Now, because we rely on RAII to release the smart pointers and destroy this memory resource on this struct,
  // the order of these entries **does** matter!
  std::pmr::monotonic_buffer_resource memory_resource;
#endif
  exprtk::parser<typename ExprTkParser::RangeField> parser{};
  exprtk::expression<typename ExprTkParser::RangeField> expression{};
  exprtk::symbol_table<typename ExprTkParser::RangeField> symbol_table{};
  std::vector<std::shared_ptr<exprtk::ifunction<typename ExprTkParser::RangeField>>> functions{};
};

ExprTkParser::ExprTkParser()
  : _data{ new ParserData{}, [](void* ptr) { delete static_cast<ParserData*>(ptr); } }
{
}

ExprTkParser::~ExprTkParser() {}

void
ExprTkParser::define_constant(const std::string& symbol, const RangeField& value)
{
  static_cast<ParserData*>(_data.get())->symbol_table.add_constant(symbol, value);
};

namespace Impl {
template<class F>
struct ExprTkFunction;

template<class R, class... Args>
struct ExprTkFunction<fu2::unique_function<R(Args...) const>>
  : public exprtk::ifunction<typename ExprTkParser::RangeField>
  , private fu2::unique_function<R(Args...) const>
{

  ExprTkFunction(std::size_t arg_size, fu2::unique_function<R(Args...) const>&& f)
    : exprtk::ifunction<typename ExprTkParser::RangeField>(arg_size)
    , fu2::unique_function<R(Args...) const>(std::move(f))
  {
    exprtk::disable_has_side_effects(*this);
  }

  inline R operator()(const Args&... args) override {
    return fu2::unique_function<R(Args...) const>::operator()(args...);
  }
};

template<class FunctionID>
void
define_function(std::size_t arg_size,
                auto& raw_data,
                const std::string& symbol,
                FunctionID&& function)
{
  ParserData& data = *static_cast<ParserData*>(raw_data.get());
  // try to allocate functions together for less jumps in memory
#ifdef DUNE_COPASI_PARSER_EXPRTK_USES_MEMORY_RESOURCE
  auto alloc = std::pmr::polymorphic_allocator<ExprTkFunction<FunctionID>>(&data.memory_resource);
#else
  auto alloc = std::allocator<ExprTkFunction<FunctionID>>();
#endif
  auto ptr = std::allocate_shared<ExprTkFunction<FunctionID>>(alloc, arg_size, std::forward<FunctionID>(function));
  data.symbol_table.add_function(symbol, *ptr);
  data.functions.emplace_back(std::move(ptr));
}

} // namespace Impl

void
ExprTkParser::define_function(const std::string& symbol, Function0D&& function)
{
  Impl::define_function(0, _data, symbol, std::move(function));
}

void
ExprTkParser::define_function(const std::string& symbol, Function1D&& function)
{
  Impl::define_function(1, _data, symbol, std::move(function));
}

void
ExprTkParser::define_function(const std::string& symbol, Function2D&& function)
{
  Impl::define_function(2, _data, symbol, std::move(function));
}

void
ExprTkParser::define_function(const std::string& symbol, Function3D&& function)
{
  Impl::define_function(3, _data, symbol, std::move(function));
}

void
ExprTkParser::define_function(const std::string& symbol, Function4D&& function)
{
  Impl::define_function(4, _data, symbol, std::move(function));
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

} // namespace Dune::Copasi
