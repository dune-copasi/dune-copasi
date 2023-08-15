#ifndef DUNE_COPASI_FILESYSTEM_HH
#define DUNE_COPASI_FILESYSTEM_HH

#include <dune-copasi-config.hh>

#if HAVE_GHC_FILESYSTEM

#include <ghc/filesystem.hpp>

namespace fs = ghc::filesystem;

#else

#include <filesystem>

namespace fs = std::filesystem;

#endif

#endif // DUNE_COPASI_FILESYSTEM_HH
