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

  dune_declare_test_label(LABELS documentation)

  file(COPY "docusaurus/docs/assets/ini"      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/docs/")

  dune_add_test(
    NAME                docs_test
    LABELS              docs
    TARGET              Dune::Copasi::Executable
    WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/docs/ini/"
    CMD_ARGS            --config=test.ini)

  if(ExprTk IN_LIST DUNE_COPASI_PARSER_LIST)
    dune_add_test(
      NAME                docs_heat
      LABELS              docs
      TARGET              Dune::Copasi::Executable
      WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/docs/ini/"
      CMD_ARGS            --config=heat.ini)
  endif()

    dune_add_test(
    NAME                docs_volka_terra
    LABELS              docs
    TARGET              Dune::Copasi::Executable
    WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/docs/ini/"
    CMD_ARGS            --config=volka_terra.ini)

  if(ExprTk IN_LIST DUNE_COPASI_PARSER_LIST)
    dune_add_test(
      NAME                docs_grey_scott
      LABELS              docs
      TARGET              Dune::Copasi::Executable
      WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/docs/ini/"
      CMD_ARGS            --config=grey_scott.ini
                          --model.time_step_operator.time_end=1.5)
  endif()

  file(COPY "docusaurus/tutorials/assets/ini" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tutorials/")

  if(ExprTk IN_LIST DUNE_COPASI_PARSER_LIST)
    dune_add_test(
      NAME                tutorial_EP_simulation1
      LABELS              docs
      TARGET              Dune::Copasi::Executable
      WORKING_DIRECTORY   "${CMAKE_CURRENT_BINARY_DIR}/tutorials/ini/"
      CMD_ARGS            --config=cardiac_ep_1.ini
                          --model.time_step_operator.time_end=10
                          --grid.refinement_level=6)
  endif()

endif()
