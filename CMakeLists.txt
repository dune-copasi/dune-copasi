cmake_minimum_required(VERSION 3.13)
project(dune-copasi CXX)

cmake_policy(SET CMP0074 NEW)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 OLD)
endif()

########## Setup dependencies ##########

# find dune dependencies
find_package(dune-common REQUIRED)

# include macros for dune projects
list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
include(DuneMacros)


# start a dune project with information from dune.module, and find dune dependencies
dune_project()

# make sure to find and install our own cmake modules
add_subdirectory(cmake/modules)

# find specific dune-copasi dependencies
find_package(TIFF REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)
find_package(Filesystem)
find_package(Eigen3)

find_package(SymEngine 0.8.0)
find_package(ExprTk)

if(NOT (SymEngine_FOUND OR ExprTk_FOUND))
  set(MUPARSER_RQUIREMENT "REQUIRED") # we need at least one parser!
endif()
find_package(muparser ${MUPARSER_RQUIREMENT})

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
    FetchContent_Declare(
      ghc_filesystem
      GIT_REPOSITORY https://github.com/gulrak/filesystem.git
      GIT_TAG        0f6f0af4b42480ad5e0e4c033e5c37659842c252
    )

    if(NOT ghc_filesystem_POPULATED)
      message("-- Populating GHC Filesystem")
      FetchContent_Populate(ghc_filesystem)
      add_subdirectory("${ghc_filesystem_SOURCE_DIR}" "${ghc_filesystem_BINARY_DIR}")
    endif()
  endif()
endif()

########## Create main targets ##########

# create main library target
dune_add_library(dune-copasi INTERFACE EXPORT_NAME Copasi)

########## Include other scripts  ##########

# include header files
add_subdirectory(dune)

# source files for final libraries and executables
add_subdirectory(src)

# documentation files
add_subdirectory(doc)

########## Setup config file ##########

# finalize the dune project: Generate and install config.h, cmake target and version file
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

# NOTE: Dowstream projects shall not generate and include config files from
# `finalize_dune_project`, otherwise, there will be repeated definitions.

# install generated config in build and install directories as <dune-copasi-config.h>
file(RENAME "${CMAKE_CURRENT_BINARY_DIR}/config.h" "${CMAKE_CURRENT_BINARY_DIR}/dune-copasi-config.h"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dune-copasi-config.h"
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

# include tests as a sub project
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# add definitions for packaging
include(SetupDuneCopasiCPack)
