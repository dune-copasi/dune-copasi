
#include <dune/copasi/common/ostream_redirect.hh>

#include <function2/function2.hpp>

#include <cstring>
#include <ios>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>

namespace Dune::Copasi {

OStreamRedirect::OStreamRedirect(std::ostream& ostream,
                                 fu2::unique_function<void(std::string_view)> redirection,
                                 bool line_buffered)
  : _ostream{ &ostream }
  , _ostream_buffer{ _ostream->rdbuf(this) } // redirect output stream buffer to this object
  , _redirection{ std::move(redirection) }
  , _line_buffered{ line_buffered }
{
  if (not _redirection) {
    throw std::ios_base::failure("stream redirection functor is empty");
  }
}

OStreamRedirect::~OStreamRedirect()
{
  try {
    // restore output stream to its original buffer
    _ostream->rdbuf(_ostream_buffer);
  } catch (const std::ios_base::failure& exception) {
    if (_redirection) {
      const auto message =
        std::string{ "OStreamRedirect was unable to restore the provided stream: " };
      _redirection(message + exception.what());
    }
  }
}

[[nodiscard]] auto
OStreamRedirect::make(std::ostream& ostream,
                      fu2::unique_function<void(std::string_view)> redirection,
                      bool line_buffered) -> std::unique_ptr<OStreamRedirect>
{
  auto* ptr = new OStreamRedirect{ ostream, std::move(redirection), line_buffered };
  return std::unique_ptr<OStreamRedirect>(ptr);
}

auto
OStreamRedirect::sync() -> int
{
  // Get a view of the currently available data
  const std::string_view buf(pbase(), pptr() - pbase());

  if (not buf.empty()) {
    // Log one message per line of output
    typename std::string_view::size_type first = 0;
    for (auto last = buf.find_first_of('\n', first); last != std::string_view::npos;
         last = buf.find_first_of('\n', first)) {
      // Ignore empty lines in unbuffered mode
      if (_line_buffered or first < last) {
        _redirection(buf.substr(first, last - first));
      }
      // skip the current character, which is a newline
      first = last + 1;
    }

    if (_line_buffered) {
      if (first == buf.size()) {
        // We have logged the entire buffer and can clear it
        seekpos(0, std::ios::out);
      } else if (first > 0 and buf.size() > 1024) {
        // There is still data in the buffer, and the buffer starts to grow a little large
        // Purge logged data from the buffer
        std::memmove(pbase(), pbase() + first, buf.size() - first);
        seekpos(buf.size() - first, std::ios::out);
      }
      // else do nothing
    } else {
      // Always print out incomplete lines as well
      if (first + 1 < buf.size() or buf[buf.size() - 1] != '\n') {
        _redirection(buf.substr(first));
      }

      // clear buffer
      seekpos(0, std::ios::out);
    }
  }
  // Forward to base class
  return Base::sync();
}

} // namespace Dune::Copasi
