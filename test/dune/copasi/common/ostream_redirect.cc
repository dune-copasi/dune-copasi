
#include <dune/copasi/common/ostream_to_spdlog.hh>

#include <spdlog/spdlog.h>

#include <iostream>

int
main(int argc, char** argv)
{

  std::cout << "Foo" << std::endl;
  {
    std::stringbuf sb;
    auto oldbuff = std::cerr.rdbuf(&sb);
    std::cerr << "stringbuf" << std::endl;
    std::cerr.rdbuf(oldbuff);
    std::cout << "Buffer contents: " << &sb;
  }
  {
    auto sr = Dune::Copasi::OStreamRedirect::make(std::cerr, [](auto sv) { spdlog::info(sv); });
    std::cerr << "raw cerr redirect" << std::endl;
  }
  {
    auto guard = Dune::Copasi::ostream2spdlog();
    std::cout << "redirect cout" << std::endl;
    std::cerr << "redirect cerr" << std::endl;
  }
}
