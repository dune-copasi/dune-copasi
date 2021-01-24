cmake_minimum_required(VERSION 3.13)
project(dune-copasi CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
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

dune_enable_all_packages()

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

target_include_directories(dune-copasi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(dune-copasi PUBLIC ${DUNE_LIBS} TIFF::TIFF muparser::muparser)

option(DUNE_COPASI_COMPILE_3D "Compile 3D cases in libraries and executables" OFF)

if (DUNE_COPASI_COMPILE_3D)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_COMPILE_3D)
endif()

if(USE_FALLBACK_FILESYSTEM)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_USE_FALLBACK_FILESYSTEM)
  target_link_libraries(dune-copasi PUBLIC ghc_filesystem)
elseif(Filesystem_FOUND)
  target_link_libraries(dune-copasi PUBLIC std::filesystem)
endif()


add_subdirectory(dune)
add_subdirectory(lib)
add_subdirectory(src)
if (dune-testtools_FOUND)
  add_subdirectory(python)
endif()
add_subdirectory(doc)
add_subdirectory(cmake/modules)

install(TARGETS dune-copasi
  EXPORT dune-copasi-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(EXPORT dune-copasi-targets
  FILE
    dune-copasi-targets.cmake
  NAMESPACE
    dune::
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/dune-copasi"
)

# finalize the dune project: Generate and install config.h, cmake targets and version
include(SetupDuneCopasiConfig)
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

add_subdirectory(test)

# include(CPack)
