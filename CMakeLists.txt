cmake_minimum_required(VERSION 3.16)
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
find_package(TBB CONFIG REQUIRED)
foreach(module ${ALL_DEPENDENCIES})
  find_package(\${module})
endforeach()
get_filename_component(_dir \"\${CMAKE_CURRENT_LIST_FILE}\" PATH)
include(\"\${_dir}/dune-copasi-installed-targets.cmake\")")

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
