#ifndef DUNE_COPASI_COMMON_OSTREAM_REDIRECT_HH
#define DUNE_COPASI_COMMON_OSTREAM_REDIRECT_HH

#include <dune-copasi-config.h>

#include <functional>
#include <memory>
#include <ostream>
#include <sstream>

namespace Dune::Copasi {

/**
 * @brief Helper class to redirect std::ostream to a functor
 * @details This object takes an `std::ostream` and redirects its input to a
 * stream buffer owned by this class. When the stream syncs, it forwards the
 * resulting strings views to a functor. It has the effect to redirect the
 * output stream to the functor.
 *
 * This `std::ostream` does not directly output data, but forwards any input to
 * its associated functor. It can operate in two different modes:
 *
 * - In line-buffered mode, the buffer will only output complete lines that have
 * been terminated by a newline character. In particular, this mode will ignore
 * any explicit flushing of the C++ stream. This mode is capable of exactly
 * reproducing the original output layout as designed by the user of the
 * `std::ostream`, but messages may appear later than expected when the user
 *   explicitly flushes the C++ stream.
 *
 * - In unbuffered mode, the buffer will always forward all pending data
 * every time the C++ stream is flushed. As most logging sinks are line-oriented
 * and insert an additional newline after each log message, this will not
 * correctly reproduce the original layout of the output. As a lot of people use
 * `std::endl` instead of just `"\n"` for ending their lines, this mode will not
 * forward empty lines to the logging system to avoid empty lines after every
 * regular line printed to the C++ stream.
 *
 * \author Steffen Müthing
 * \copyright BSD-2-Clause
 */
class OStreamRedirect final : public std::stringbuf
{

  using Base = std::stringbuf;

  //! internal constructor of the redirection object
  OStreamRedirect(std::ostream& ostream,
                  std::move_only_function<void(std::string_view)> redirection,
                  bool line_buffered);

public:
  ~OStreamRedirect() override;

  /**
   * @brief RAII safe constructor of the redirection object
   * While the object is alive, the output of ostream will be redirected to the
   * redirection functor
   *
   * @param ostream        The output stream to redirect
   * @param redirection    The functor to redirect the stream of characters
   * @param line_buffered  Whether to operate on buffered mode
   * @return RAII safe pointer to a redirection object
   */
  [[nodiscard]] static std::unique_ptr<OStreamRedirect> make(
    std::ostream& ostream,
    std::move_only_function<void(std::string_view)> redirection,
    bool line_buffered = true);

  OStreamRedirect(const OStreamRedirect& other) = delete;
  OStreamRedirect& operator=(const OStreamRedirect& other) = delete;

  OStreamRedirect(OStreamRedirect&& other) noexcept = delete;
  OStreamRedirect& operator=(OStreamRedirect&& other) noexcept = delete;

private:
  int sync() override;
  std::ostream* _ostream;
  std::streambuf* _ostream_buffer;
  std::move_only_function<void(std::string_view)> _redirection;
  bool _line_buffered;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_OSTREAM_REDIRECT_HH
