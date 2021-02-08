cmake_minimum_required(VERSION 3.13)
project(dune-copasi CXX)

include(GNUInstallDirs)

cmake_policy(SET CMP0074 NEW)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 OLD)
endif()

# define cmake options
option(DUNE_COPASI_SD_EXECUTABLE
       "Build SingleDomain executable by default"
       OFF
)
option(DUNE_COPASI_SD_LIBRARY
       "Include precompiled SingleDomain objects in dune-copasi library (optional)"
       ${DUNE_COPASI_SD_EXECUTABLE}
)
option(DUNE_COPASI_MD_EXECUTABLE
       "Build MultiDimain executable by default"
       ON
)
option(DUNE_COPASI_MD_LIBRARY
       "Include precompiled MultiDimain objects in dune-copasi library (optional)"
       ${DUNE_COPASI_MD_EXECUTABLE}
)
option(DUNE_COPASI_COMPILE_3D
      "Compile 3D cases in libraries and executables"
      OFF
)

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
find_package(muparser REQUIRED)
find_package(TIFF REQUIRED)
find_package(Filesystem)

# if standard C++ filesystem is not found, using a fallback is mandatory
include(CMakeDependentOption)
cmake_dependent_option(DUNE_USE_FALLBACK_FILESYSTEM
    "Use bundled std::filesystem alternative (enable this for macOS older than 10.15)" OFF
    "Filesystem_FOUND" ON
)

if(DUNE_USE_FALLBACK_FILESYSTEM)
  # let's use ghc_filesystem
  find_package(ghc_filesystem)

  # ... library not installed, we provide it from git
  if (NOT ghc_filesystem_FOUND)
    include(FetchContent)
    message("-- Declaring GHC Filesystem")
    FetchContent_Declare(
      ghc_filesystem
      GIT_REPOSITORY https://github.com/gulrak/filesystem
      GIT_TAG        v1.5.0
    )

    if(NOT ghc_filesystem_POPULATED)
      message("-- Populating GHC Filesystem")
      FetchContent_Populate(ghc_filesystem)
      add_subdirectory("${ghc_filesystem_SOURCE_DIR}" "${ghc_filesystem_BINARY_DIR}")
      export(EXPORT ghcFilesystemConfig FILE "${CMAKE_CURRENT_BINARY_DIR}/ghc_filesystem-targets.cmake")
    endif()
  endif()
endif()

# create main target
add_library(dune-copasi)
set_target_properties(dune-copasi PROPERTIES LIBRARY_OUTPUT_NAME dunecopasi)
set_target_properties(dune-copasi PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasi)

# include header files
add_subdirectory(dune)

# include main library set up
add_subdirectory(lib)

# source files for final libraries and executables
add_subdirectory(src)

# documentation files
add_subdirectory(doc)

# finalize the dune project: Generate and install config.h, cmake target and version file
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
target_compile_definitions(dune-copasi INTERFACE HAVE_DUNE_COPASI_CONFIG_H)

# main library is completely set up, define a read only alias for our main target
add_library(dune-copasi::dune-copasi ALIAS dune-copasi)

# install generated config in build and install directories
file(COPY "${CMAKE_CURRENT_BINARY_DIR}/config.h"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dune/copasi"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dune/copasi/config.h"
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi"
)

# If the FC.h file is generated (dune-common<=2.7), install into build and install directories
file(TOUCH "${CMAKE_CURRENT_BINARY_DIR}/FC.h")
file(COPY "${CMAKE_CURRENT_BINARY_DIR}/FC.h"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dune/copasi"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dune/copasi/FC.h"
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi"
)

# taget configuration for downstream projects
install(TARGETS dune-copasi
        EXPORT dune-copasi-targets
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          COMPONENT Libraries
          NAMELINK_COMPONENT Development
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
          COMPONENT Development
)

# export targets to be used with installed tree
install(EXPORT dune-copasi-targets
        FILE dune-copasi-targets.cmake
        NAMESPACE dune-copasi::
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dune-copasi"
)

# export targets to be used with build tree
export(EXPORT dune-copasi-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/dune-copasi-targets.cmake"
       NAMESPACE dune-copasi::
)

# include tests as a sub project
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# add definitions for packaging
include(SetupDuneCopasiCPack)