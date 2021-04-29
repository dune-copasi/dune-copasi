if(DUNE_COPASI_INSTALL_XEUS_CLING)

  # get module that makes the set up
  include(FetchContent)
  message("-- Declaring Xeus Cling Setup")
  FetchContent_Declare(
    xeus_cling_setup
    GIT_REPOSITORY https://github.com/dokempf/xeus-cling-cmake-setup.git
    GIT_TAG        6356d6ce997f29726463d1422ed477c0b4539396
  )

  if(NOT xeus_cling_setup_POPULATED)
    message("-- Populating Xeus Cling Setup")
    FetchContent_Populate(xeus_cling_setup)
    add_subdirectory("${xeus_cling_setup_SOURCE_DIR}" "${xeus_cling_setup_BINARY_DIR}")
  endif()

  if(TARGET dune-logging-fmt)
    get_target_property(xeus_cling_include dune-logging-fmt INTERFACE_INCLUDE_DIRECTORIES)
  endif()

  xeus_cling_setup(
    TARGETS dune-copasi::dune-copasi
    INCLUDE_DIRECTORIES "${xeus_cling_include}"
    COMPILE_FLAGS
      -DENABLE_UG=1
      -DUG_USE_NEW_DIMENSION_DEFINES
    SETUP_HEADERS
      ${CMAKE_BINARY_DIR}/dune/copasi/config.h
    KERNEL_NAME "C++17 (DUNE-COPASI)"
    REQUIRED
  )

  # install kernel by default
  add_dependencies(install_kernelspec dune-copasi)
  set_target_properties(install_kernelspec PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()
