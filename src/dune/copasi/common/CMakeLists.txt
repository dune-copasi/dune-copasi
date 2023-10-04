
# find specific dune-copasi dependencies
find_package(TIFF REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)
find_package(Filesystem)

# if standard C++ filesystem is not found, using a fallback is mandatory
include(CMakeDependentOption)
cmake_dependent_option(DUNE_USE_FALLBACK_FILESYSTEM
    "Use bundled std::filesystem alternative (enable this for macOS older than 10.15)" OFF
    "Filesystem_FOUND" ON
)

if(DUNE_USE_FALLBACK_FILESYSTEM)
  # let's use ghc_filesystem
  find_package(ghc_filesystem QUIET)

  # ... library not installed, we provide it from git
  if (NOT ghc_filesystem_FOUND)
    include(FetchContent)
    message("-- Declaring GHC Filesystem")
    option(GHC_FILESYSTEM_WITH_INSTALL "With install target" ON)
    FetchContent_Declare(
      ghc_filesystem
      GIT_REPOSITORY https://github.com/gulrak/filesystem.git
      GIT_TAG        v1.5.14
    )

    if(NOT ghc_filesystem_POPULATED)
      message("-- Populating GHC Filesystem")
      FetchContent_Populate(ghc_filesystem)
      add_subdirectory("${ghc_filesystem_SOURCE_DIR}" "${ghc_filesystem_BINARY_DIR}")
    endif()
  endif()
endif()

# The function2 library provides an alternative to std::move_only_function
# (which otherwise requires c++23 in libstdc++ and is not yet implemented in libc++)
include(FetchContent)
message("-- Declaring function2")
FetchContent_Declare(
        function2
        GIT_REPOSITORY https://github.com/Naios/function2
        GIT_TAG        4.2.3
)
if(NOT function2_POPULATED)
  message("-- Populating function2")
  FetchContent_Populate(function2)
  add_library(function2 INTERFACE)
  install(DIRECTORY ${function2_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/external/function2)
  target_include_directories(function2
    INTERFACE
      $<BUILD_INTERFACE:${function2_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/dune/external/function2>)
  install(TARGETS function2 EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
endif()

# include dependencies in config file
set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_package(TIFF REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)
")

if(${DUNE_USE_FALLBACK_FILESYSTEM})
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}\nfind_package(ghc_filesystem REQUIRED)")
else()
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}\nfind_package(Filesystem REQUIRED)")
endif()

set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}" PARENT_SCOPE)

add_library(dune-copasi-common STATIC)
add_library(Dune::Copasi::Common ALIAS dune-copasi-common)
set_target_properties(dune-copasi-common PROPERTIES EXPORT_NAME Copasi::Common)

target_sources(dune-copasi-common
  PRIVATE
    tiff_file.cc
    tiff_grayscale.cc
    ostream_redirect.cc
)

target_include_directories(dune-copasi-common
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${PROJECT_BINARY_INCLUDEDIR}>
)

target_compile_definitions(dune-copasi-common
  PUBLIC
    $<$<BOOL:${DUNE_USE_FALLBACK_FILESYSTEM}>:HAVE_GHC_FILESYSTEM>
)

target_link_libraries(dune-copasi-common
  PUBLIC
    $<TARGET_NAME_IF_EXISTS:ghcFilesystem::ghc_filesystem>
    function2
  PRIVATE
    TIFF::TIFF
  INTERFACE
    fmt::fmt
    spdlog::spdlog
    ${DUNE_LIBS}
)

target_compile_features(dune-copasi-common PUBLIC cxx_std_20)

# target configuration for downstream projects
install(TARGETS dune-copasi-common EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
