cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(dune-copasi CXX)

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

########## Include other scripts  ##########

# include header files
add_subdirectory(dune)

# source files for final libraries and executables
add_subdirectory(src)

# documentation files
add_subdirectory(doc)

########## Setup config file ##########

# include targets file in config file
set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
set(DUNE_COPASI_GRID_DIMENSIONS ${DUNE_COPASI_GRID_DIMENSIONS})
set(DUNE_COPASI_PARSER_LIST ${DUNE_COPASI_PARSER_LIST})
# ensure dune cmake modules are found and used
set(PREVIOUS_CMAKE_MODULE_PATH \${CMAKE_MODULE_PATH})
list(PREPEND CMAKE_MODULE_PATH \"\${dune-copasi_MODULE_PATH}\")

include(CMakeFindDependencyMacro)
")

foreach(module ${DUNE_FOUND_DEPENDENCIES})
  if(${${module}_FOUND})
    set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(${module} REQUIRED)
list(APPEND CMAKE_MODULE_PATH \"\${${module}_MODULE_PATH}\")")
  endif()
endforeach()


if(TARGET MPI::MPI_C)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(MPI 3.0 REQUIRED COMPONENTS C)")
endif()

set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(TIFF REQUIRED)
find_dependency(spdlog REQUIRED)
find_dependency(fmt REQUIRED)
find_dependency(TBB CONFIG REQUIRED)")

if(TARGET METIS::METIS)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(METIS REQUIRED)")
endif()

if(TARGET Dune::Wrapper::SymEngine)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(SymEngine 0.8.0 REQUIRED)
add_library(Dune::Wrapper::SymEngine INTERFACE IMPORTED)
target_link_libraries(Dune::Wrapper::SymEngine INTERFACE \${SYMENGINE_LIBRARIES})
target_include_directories(Dune::Wrapper::SymEngine INTERFACE \${SYMENGINE_INCLUDE_DIRS})")
endif()

if(TARGET muparser::muparser)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(muparser REQUIRED)")
endif()

if(TARGET parafields::parafields)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(parafields REQUIRED)")
endif()

if(TARGET SuiteSparse::UMFPACK)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(SuiteSparse REQUIRED COMPONENTS UMFPACK)")
endif()

if(TARGET SuperLU::SuperLU)
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}
find_dependency(SuperLU REQUIRED)")
endif()

set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}

list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
get_filename_component(_dir \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)
include(\"\${_dir}/dune-copasi-installed-targets.cmake\")
# restore previous CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH \${PREVIOUS_CMAKE_MODULE_PATH})
unset(PREVIOUS_CMAKE_MODULE_PATH)
")

# finalize the dune project: Generate and install config headers, cmake targets and version files
finalize_dune_project()

# export targets to be used with installed tree
install(EXPORT dune-copasi-targets
        FILE dune-copasi-installed-targets.cmake
        NAMESPACE Dune::
        COMPONENT Development
        DESTINATION "${DUNE_INSTALL_LIBDIR}/cmake/dune-copasi"
)

# export targets to be used with build tree
export(EXPORT dune-copasi-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/dune-copasi-installed-targets.cmake"
       NAMESPACE Dune::
)
# include tests as a sub project
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# add definitions for packaging
include(SetupDuneCopasiCPack)
