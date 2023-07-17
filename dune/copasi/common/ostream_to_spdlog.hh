#ifndef DUNE_COPASI_COMMON_OSTREAM_TO_SPDLOG_HH
#define DUNE_COPASI_COMMON_OSTREAM_TO_SPDLOG_HH

#include <dune/copasi/common/ostream_redirect.hh>

#include <spdlog/spdlog.h>

#include <array>
#include <iostream>

namespace Dune::Copasi {

/**
 * @brief Standard output redirection to spdlog
 * @details During the lifetime of the obtained object, it will redirect every
 * character streamed to std::cout, std::cerr, and std::clog towards the
 * spdlog logging system
 */
[[nodiscard]] inline auto
ostream2spdlog()
{
  return std::array{ OStreamRedirect::make(std::cout, [](auto msg) { spdlog::info(msg); }),
                     OStreamRedirect::make(std::cerr, [](auto msg) { spdlog::critical(msg); }),
                     OStreamRedirect::make(std::clog, [](auto msg) { spdlog::debug(msg); }) };
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_OSTREAM_TO_SPDLOG_HH
