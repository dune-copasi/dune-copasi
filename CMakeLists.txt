cmake_minimum_required(VERSION 3.13)
project(dune-copasi CXX)

cmake_policy(SET CMP0074 NEW)

# find dune-common
find_package(dune-common REQUIRED)

# make sure to find cmake modules
list(APPEND CMAKE_MODULE_PATH
  "${PROJECT_SOURCE_DIR}/cmake/modules"
  "${dune-common_MODULE_PATH}"
)

# include macros for dune projects
include(DuneMacros)

option(USE_FALLBACK_FILESYSTEM
      "Use bundled std::filesystem alternative (enable this for macOS older than 10.15)"
      OFF
)

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

# find dependencies
find_package(muparser REQUIRED)
find_package(TIFF REQUIRED)
find_package(Filesystem)
if(USE_FALLBACK_FILESYSTEM)
  find_package(ghc_filesystem)
endif()

# start a dune project with information from dune.module
dune_project()

# main target
add_library(dune-copasi)
set_target_properties(dune-copasi PROPERTIES LIBRARY_OUTPUT_NAME dunecopasi)
set_target_properties(dune-copasi PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasi)

# add dune dependencies to this target
dune_target_enable_all_packages(dune-copasi)

# add dependencies
target_link_libraries(dune-copasi PUBLIC ${DUNE_LIBS} TIFF::TIFF muparser::muparser)

if (DUNE_COPASI_COMPILE_3D)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_COMPILE_3D)
endif()

add_subdirectory(dune)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(doc)

# finalize the dune project: Generate and install config.h, cmake target and version file
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

# set a read only alias for our main target
add_library(dune::copasi ALIAS dune-copasi)

# include tests as a sub project
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# add definitions for packaging
include(SetupDuneCopasiCPack)