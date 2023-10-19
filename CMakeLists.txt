cmake_minimum_required(VERSION 3.13)
project(dune-copasi CXX)

include(GNUInstallDirs)

cmake_policy(SET CMP0074 NEW)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 OLD)
endif()

# define cmake options
option(DUNE_COPASI_SD_EXECUTABLE
       "Build and install SingleDomain executable"
       OFF
)
option(DUNE_COPASI_MD_EXECUTABLE
       "Build and install MultiDimain executable"
       ON
)
option(DUNE_COPASI_COMPILE_3D
      "Compile 3D cases in executables"
      OFF
)
option(DUNE_COPASI_INSTALL_XEUS_CLING
      "Setup just in time compilation to enable jupyter notebooks (only with build tree)"
      OFF
)
option(DUNE_COPASI_USE_MP_UNITS
      "Use types to represent units with the mp-units library"
      OFF
)
########## Setup dependencies ##########

# find dune dependencies
find_package(dune-common REQUIRED)

# include macros for dune projects
list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
include(DuneMacros)


# start a dune project with information from dune.module, and find dune dependencies
dune_project()

# find_package(dune-assembler REQUIRED)

# make sure to find and install our own cmake modules
add_subdirectory(cmake/modules)

# find specific dune-copasi dependencies
find_package(TIFF REQUIRED)
find_package(Filesystem)

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

if (DUNE_COPASI_USE_MP_UNITS)
find_package(gsl-lite CONFIG REQUIRED)
  if (gsl-lite_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      mp_units
      GIT_REPOSITORY https://github.com/mpusz/units.git
      GIT_TAG        v0.7.0
    )

    if(NOT mp_units_POPULATED)
      FetchContent_Populate(mp_units)
      add_subdirectory("${mp_units_SOURCE_DIR}/src" "${mp_units_SOURCE_DIR}/build")
    endif()
  endif()
endif()

########## Create main targets ##########

# create main library target
add_library(dune-copasi)
set_property(TARGET dune-copasi PROPERTY EXPORT_NAME Copasi)
add_library(Dune::Copasi ALIAS dune-copasi)
set_target_properties(dune-copasi PROPERTIES LIBRARY_OUTPUT_NAME dunecopasi)
set_target_properties(dune-copasi PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasi)

# add_dune_alberta_flags(dune-copasi WORLDDIM 2)
add_dune_alberta_flags(dune-copasi WORLDDIM 3)

# get properties defined by dune macros
get_property(dune_all_libs GLOBAL PROPERTY ALL_PKG_LIBS)
get_property(dune_all_incs GLOBAL PROPERTY ALL_PKG_INCS)
get_property(dune_all_defs GLOBAL PROPERTY ALL_PKG_DEFS)
get_property(dune_all_opts GLOBAL PROPERTY ALL_PKG_OPTS)


target_link_libraries(dune-copasi
  PUBLIC
    ${DUNE_LIBS}
    Dune::Assembler
    TIFF::TIFF
  INTERFACE
    ${dune_all_libs}
    $<$<BOOL:${muparser_FOUND}>:muparser::muparser>
    $<$<BOOL:${SymEngine_FOUND}>:${SYMENGINE_LIBRARIES}>
    $<$<BOOL:${ExprTk_FOUND}>:ExprTk::ExprTk>
    $<IF:$<BOOL:${DUNE_USE_FALLBACK_FILESYSTEM}>,ghcFilesystem::ghc_filesystem,std::filesystem>
    $<$<BOOL:${DUNE_COPASI_USE_MP_UNITS}>:mp-units::mp-units>
)

target_include_directories(dune-copasi
  INTERFACE
    ${dune_all_incs}
    $<$<BOOL:${SymEngine_FOUND}>:${SYMENGINE_INCLUDE_DIRS}>)

target_compile_definitions(dune-copasi
  INTERFACE
    ${dune_all_defs}
    HAVE_MUPARSER=$<BOOL:${muparser_FOUND}>
    HAVE_SYMENGINE=$<BOOL:${SymEngine_FOUND}>
    HAVE_EXPRTK=$<BOOL:${ExprTk_FOUND}>
    HAVE_GHC_FILESYSTEM=$<BOOL:${DUNE_USE_FALLBACK_FILESYSTEM}>
    HAVE_UNITS=$<BOOL:${DUNE_COPASI_USE_MP_UNITS}>
)

target_compile_options(dune-copasi
  INTERFACE
    ${dune_all_opts})

# create executable targets
add_executable(singledomain-exec)
add_executable(Dune::Copasi::SingleDomainExec ALIAS singledomain-exec)

add_executable(multidomain-exec)
add_executable(Dune::Copasi::MultiDomainExec ALIAS multidomain-exec)

add_executable(multidomain-membrane-exec)
add_executable(Dune::Copasi::MultiDomainMembraneExec ALIAS multidomain-membrane-exec)

########## Include other scripts  ##########

# include header files
add_subdirectory(dune)

# source files for final libraries and executables
add_subdirectory(src)

# jupyter notebooks
# add_subdirectory(notebooks)

# documentation files
add_subdirectory(doc)

########## Setup config file ##########

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

# taget configuration for downstream projects
install(TARGETS dune-copasi
        EXPORT dune-copasi-targets
        LIBRARY DESTINATION "${DUNE_INSTALL_LIBDIR}"
          COMPONENT Libraries
          NAMELINK_COMPONENT Development
        ARCHIVE DESTINATION "${DUNE_INSTALL_LIBDIR}"
          COMPONENT Development
)

# export targets to be used with installed tree
install(EXPORT dune-copasi-targets
        FILE dune-copasi-targets.cmake
        NAMESPACE Dune::
        COMPONENT Development
        DESTINATION "${DUNE_INSTALL_LIBDIR}/cmake/dune-copasi"
)

# export targets to be used with build tree
export(EXPORT dune-copasi-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/dune-copasi-targets.cmake"
       NAMESPACE Dune::
)

# include tests as a sub project
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# add definitions for packaging
include(SetupDuneCopasiCPack)
