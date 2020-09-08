# We require version CMake version 3.1 to prevent issues
# with dune_enable_all_packages and older CMake versions.
cmake_minimum_required(VERSION 3.1)
project(dune-copasi CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

cmake_policy(SET CMP0074 NEW)

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages()

set(USE_FALLBACK_FILESYSTEM no CACHE BOOL "Use bundled std::filesystem alternative (enable this for macOS older than 10.15)")
if(USE_FALLBACK_FILESYSTEM)
  add_subdirectory(ext/ghc)
  install(FILES ext/ghc/include/ghc/filesystem.hpp
          COMPONENT Development
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ghc")
endif()
add_subdirectory(src)
add_subdirectory(dune)
add_subdirectory(lib)
if (dune-testtools_FOUND)
  add_subdirectory(python)
endif()
add_subdirectory(doc)
add_subdirectory(test)
add_subdirectory(cmake/modules)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)