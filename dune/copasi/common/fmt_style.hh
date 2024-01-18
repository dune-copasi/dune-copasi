#ifndef DUNE_COPASI_COMMON_FMT_STYLE_HH
#define DUNE_COPASI_COMMON_FMT_STYLE_HH

#include <dune-copasi-config.hh>

#include <fmt/color.h>

#if FMT_VERSION >= 90000
#define DUNE_COPASI_FMT_STYLED_BOLD(key) fmt::styled(key, fmt::emphasis::bold)
#define DUNE_COPASI_FMT_STYLED_ITALIC(key) fmt::styled(key, fmt::emphasis::italic)
#define DUNE_COPASI_FMT_STYLED_DARK_GRAY(key) fmt::styled(key, fmt::fg(fmt::color::dark_gray))
#else
#define DUNE_COPASI_FMT_STYLED_BOLD(key) key
#define DUNE_COPASI_FMT_STYLED_ITALIC(key) key
#define DUNE_COPASI_FMT_STYLED_DARK_GRAY(key) key
#endif

#endif // DUNE_COPASI_COMMON_FMT_STYLE_HH
