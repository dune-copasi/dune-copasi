#ifndef DUNE_COPASI_FILESYSTEM_HH
#define DUNE_COPASI_FILESYSTEM_HH

#ifdef DUNE_USE_FALLBACK_FILESYSTEM

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
