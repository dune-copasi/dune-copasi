cmake_minimum_required(VERSION 3.13)
project(dune-copasi-tests CXX)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-copasi as if this were downstream project
  find_package(dune-copasi REQUIRED IMPORTED COMPONENTS Development Runtime)
endif()

add_subdirectory(dune)

enable_testing()

if (TARGET Dune::Copasi::Executable)
  foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
    dune_add_test(
      NAME            gauss-${dim}D
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/gauss.ini --grid.dimension=${dim})

    dune_add_test(
      NAME            exp-${dim}D
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/exp.ini --grid.dimension=${dim})
  endforeach()
endif()


if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  set_target_properties(build_unit_tests build_system_tests PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()
