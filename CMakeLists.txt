cmake_minimum_required(VERSION 3.13)
project(dune-copasi CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      "${PROJECT_BINARY_DIR}")
endif()

cmake_policy(SET CMP0074 NEW)

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
  "${PROJECT_SOURCE_DIR}/cmake/modules"
  "${dune-common_MODULE_PATH}"
)

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

set(USE_FALLBACK_FILESYSTEM no CACHE BOOL "Use bundled std::filesystem alternative (enable this for macOS older than 10.15)")
if(USE_FALLBACK_FILESYSTEM)
  add_subdirectory(ext/ghc)
  install(FILES ext/ghc/include/ghc/filesystem.hpp
          COMPONENT Development
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ghc")
endif()


# Single domain target
set(DUNE_COPASI_SD_EXECUTABLE OFF
    CACHE BOOL "Build SingleDomain executable by default"
)
set(DUNE_COPASI_SD_LIBRARY ${DUNE_COPASI_SD_EXECUTABLE}
    CACHE BOOL "Include precompiled SingleDomain objects in dune-copasi library (optional)"
)

# Multiple domain target
set(DUNE_COPASI_MD_EXECUTABLE ON
    CACHE BOOL "Build MultiDimain executable by default"
)
set(DUNE_COPASI_MD_LIBRARY ${DUNE_COPASI_MD_EXECUTABLE}
    CACHE BOOL "Include precompiled MultiDimain objects in dune-copasi library (optional)"
)

# main target
add_library(dune-copasi)
set_target_properties(dune-copasi PROPERTIES LIBRARY_OUTPUT_NAME dunecopasi)
set_target_properties(dune-copasi PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasi)

# add dune dependencies to this target
dune_target_enable_all_packages(dune-copasi)

# add directories for includes and libs
target_link_directories(dune-copasi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>
)

target_include_directories(dune-copasi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/lib>
  $<INSTALL_INTERFACE:include>
)

# add dependencies
target_link_libraries(dune-copasi PUBLIC ${DUNE_LIBS} TIFF::TIFF muparser::muparser)

if(USE_FALLBACK_FILESYSTEM)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_USE_FALLBACK_FILESYSTEM)
  target_link_libraries(dune-copasi PUBLIC ghc_filesystem)
elseif(Filesystem_FOUND)
  target_link_libraries(dune-copasi PUBLIC std::filesystem)
endif()

option(DUNE_COPASI_COMPILE_3D "Compile 3D cases in libraries and executables" OFF)

if (DUNE_COPASI_COMPILE_3D)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_COMPILE_3D)
endif()

add_subdirectory(dune)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

# if (DUNE_PYTHON_VIRTUALENV_SETUP)
#   dune_python_install_package(PATH python)
#   add_subdirectory(python)
# endif()

# finalize the dune project: Generate and install config.h, cmake target and version file
include(SetupDuneCopasiConfig)
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

target_compile_definitions(dune-copasi INTERFACE HAVE_DUNE_COPASI_CONFIG_H)

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

# export targets to be used with the installed package
install(EXPORT dune-copasi-targets
        FILE dune-copasi-targets.cmake
        NAMESPACE dune-copasi::
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dune-copasi"
)

# export targets to be used with build directories
export(EXPORT dune-copasi-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/dune-copasi-targets.cmake"
       NAMESPACE dune-copasi::
)

# include tests as a sub project
add_library(dune-copasi::dune-copasi ALIAS dune-copasi)
add_subdirectory(test)
