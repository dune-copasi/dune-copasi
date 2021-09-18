cmake_minimum_required(VERSION 3.9)
project(dune-copasi-notebooks)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-copasi as if this were downstream project
  find_package(dune-copasi REQUIRED IMPORTED COMPONENTS Libraries Development)
  set(DUNE_COPASI_INSTALL_XEUS_CLING TRUE)
endif()

if(DUNE_COPASI_INSTALL_XEUS_CLING)
  # get module that makes the set up
  include(FetchContent)
  message("-- Declaring Xeus Cling Setup")
  FetchContent_Declare(
    xeus_cling_setup
    GIT_REPOSITORY https://github.com/dokempf/xeus-cling-cmake-setup.git
    GIT_TAG        a8f955cc3284fbf9f66cca5e4f3cdac54d51ffa8
  )

  if(NOT xeus_cling_setup_POPULATED)
    message("-- Populating Xeus Cling Setup")
    FetchContent_Populate(xeus_cling_setup)
    add_subdirectory("${xeus_cling_setup_SOURCE_DIR}" "${xeus_cling_setup_BINARY_DIR}")
  endif()

  xeus_cling_setup(
    TARGETS Dune::Copasi
    SETUP_HEADERS
      dune/copasi/config.h
    LINK_LIBRARIES
      /opt/conda/lib/libmuparser.so # Help! need this because transitive dependencies do not work somehow
    KERNEL_NAME "C++17 (DUNE-COPASI)"
    KERNEL_LOGO_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../doc/docusaurus/static/img/logo-64x64.png"
    REQUIRED
  )

  # install kernel by default
  add_dependencies(install_kernelspec Dune::Copasi)
  set_target_properties(install_kernelspec PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()
