#ifndef DUNE_COPASI_FILESYSTEM_HH
#define DUNE_COPASI_FILESYSTEM_HH

#if HAVE_GHC_FILESYSTEM

#include <ghc/filesystem.hpp>
namespace fs {
using namespace ghc::filesystem;
}

#else

#include <filesystem>
namespace fs {
using namespace std::filesystem;
}

#endif

#endif // DUNE_COPASI_FILESYSTEM_HH
