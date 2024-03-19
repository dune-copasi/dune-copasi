cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(dune-copasi-docs CXX C)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-copasi as if this were downstream project
  find_package(dune-copasi REQUIRED IMPORTED COMPONENTS Development Runtime)
  set(DUNE_BUILD_TESTS_ON_MAKE_ALL ON)
else()
  set(DUNE_BUILD_TESTS_ON_MAKE_ALL OFF)
endif()

add_subdirectory(doxygen)
add_subdirectory(man)

if (TARGET Dune::Copasi::Executable)
  list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
  list(APPEND CMAKE_MODULE_PATH "${dune-pdelab_MODULE_PATH}")
  include(DuneMacros)
  include(DuneAddGTest)

  enable_testing()

  dune_declare_test_label(LABELS tutorial)

  file(COPY "docusaurus/tutorials/assets/ini" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tutorials/")

  dune_add_test(
    NAME                tutorial_test
    LABELS              tutorial
    TARGET              Dune::Copasi::Executable
    WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/tutorials/ini/"
    CMD_ARGS            --config=test.ini)

  if(ExprTk IN_LIST DUNE_COPASI_PARSER_LIST)
    dune_add_test(
      NAME                tutorial_EP_simulation1
      LABELS              tutorial
      TARGET              Dune::Copasi::Executable
      WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/tutorials/ini/"
      CMD_ARGS            --config=EP_simulation1.ini
                          --model.time_step_operator.time_end=10)
  endif()

endif()
