cmake_minimum_required(VERSION 3.13)
project(dune-copasi-tests CXX C)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-copasi as if this were downstream project
  find_package(dune-copasi REQUIRED IMPORTED COMPONENTS Development Runtime)
  set(DUNE_BUILD_TESTS_ON_MAKE_ALL ON)
else()
  set(DUNE_BUILD_TESTS_ON_MAKE_ALL OFF)
endif()

list(APPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
include(DuneMacros)

# add_subdirectory(dune)

enable_testing()

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data COPY_ON_ERROR SYMBOLIC)

dune_declare_test_label(LABELS unit system)

if (TARGET Dune::Copasi::Executable)

  foreach(simple_opt help;help-full;version;parser-list;dump-config)
    dune_add_test(
      NAME            ${simple_opt}
      TARGET          Dune::Copasi::Executable
      LABELS          system
      CMD_ARGS        --${simple_opt})
  endforeach()

  dune_add_test(
    NAME            empty
    TARGET          Dune::Copasi::Executable
    LABELS          system
    EXPECT_FAIL)

  if(ExprTk IN_LIST DUNE_COPASI_PARSER_LIST)
    dune_add_test(
      NAME            disk
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/disk.ini
      LABELS          system)
  endif()

  foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
    dune_add_test(
      NAME            exp-${dim}D
      TARGET          Dune::Copasi::Executable
      LABELS          system
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/exp.ini
                      --grid.dimension=${dim})

    foreach(parser ${DUNE_COPASI_PARSER_LIST})
      dune_add_test(
        NAME            gauss-${dim}D-${parser}
        TARGET          Dune::Copasi::Executable
        LABELS          system
        CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/gauss.ini
                        --grid.dimension=${dim}
                        --model.parser_type=${parser})
    endforeach()
  endforeach()
endif()
