if(DUNE_COPASI_INSTALL_XEUS_CLING)

  # get module that makes the set up
  include(FetchContent)
  message("-- Declaring Xeus Cling Setup")
  FetchContent_Declare(
    xeus_cling_setup
    GIT_REPOSITORY https://github.com/dokempf/xeus-cling-cmake-setup.git
    GIT_TAG        3f94836ea521b76dd89f24752ad12e320806f4b2
  )

  if(NOT xeus_cling_setup_POPULATED)
    message("-- Populating Xeus Cling Setup")
    FetchContent_Populate(xeus_cling_setup)
    add_subdirectory("${xeus_cling_setup_SOURCE_DIR}" "${xeus_cling_setup_BINARY_DIR}")
  endif()

  xeus_cling_setup(
    TARGETS dune-copasi::dune-copasi
    COMPILE_DEFINITIONS
      UG_USE_NEW_DIMENSION_DEFINES # TODO: remove on dune 2.8
    SETUP_HEADERS
      ${CMAKE_BINARY_DIR}/dune/copasi/config.h
    KERNEL_NAME "C++17 (DUNE-COPASI)"
    KERNEL_LOGO_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../doc/docusaurus/static/img/logo-64x64.png"
    REQUIRED
  )

  # install kernel by default
  add_dependencies(install_kernelspec dune-copasi)
  set_target_properties(install_kernelspec PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()
