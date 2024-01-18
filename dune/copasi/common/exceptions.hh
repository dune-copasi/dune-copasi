#ifndef DUNE_COPASI_COMMON_EXCEPTION_HH
#define DUNE_COPASI_COMMON_EXCEPTION_HH

#include <dune-copasi-config.hh>

#include <dune/common/classname.hh>

#include <fmt/core.h>

#include <utility>
#if __has_include(<stacktrace>)
#include <stacktrace>
#endif

#if __has_include(<format>)
#include <format>
#endif

namespace Dune::Copasi {

template<typename Exception, typename... Args>
[[nodiscard]] inline auto
format_exception(Exception&& e, fmt::format_string<Args...> format, Args&&... args)
{
  auto message = fmt::format(std::move(format), std::forward<Args>(args)...);
  message += fmt::format("\n\nException type: {}", className<Exception>());
#if ___cpp_lib_stacktrace >= 202011L && __cpp_lib_formatters >= 202302L
  message += std::format("\nStacktrace:\n{}", std::stacktrace::current());
#endif
  e.message(message);
  return std::forward<Exception>(e);
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_EXCEPTION_HH
