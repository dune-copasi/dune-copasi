cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
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
list(APPEND CMAKE_MODULE_PATH "${dune-pdelab_MODULE_PATH}")
include(DuneMacros)
include(DuneAddGTest)

add_subdirectory(dune)

enable_testing()

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data COPY_ON_ERROR SYMBOLIC)

dune_declare_test_label(LABELS unit system)

if (TARGET Dune::Copasi::Executable)

  foreach(simple_opt help;help-full;version;parser-list;dump-config)
    dune_add_test(
      NAME            ${simple_opt}
      LABELS          system
      TIMEOUT         10
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --${simple_opt})
  endforeach()

  dune_add_test(
    NAME            empty
    TARGET          Dune::Copasi::Executable
    LABELS          system
    TIMEOUT         10
    EXPECT_FAIL)

  # over|exact|under the proper time-step to snap time to the end time
  foreach(init_step 0.199999 0.100000 0.999999)
    dune_add_test(
      NAME            time-snap-${init_step}
      TIMEOUT         10
      LABELS          system
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/time_snap.ini
                      --model.time_step_operator.time_step_initial=${init_step})
  endforeach()

  if(ExprTk IN_LIST DUNE_COPASI_PARSER_LIST)
    dune_add_test(
      NAME            two-disks
      TIMEOUT         10
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/two_disks.ini
      LABELS          system)

    dune_add_test(
      NAME            two-disks-cell-data
      TIMEOUT         10
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/two_disks_cell_data.ini
      LABELS          system)
  
    dune_add_test(
      NAME            mitchell-schaefer
      LABELS          system
      TIMEOUT         10
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/mitchell_schaefer.ini
                      --grid.refinement_level=1
                      --model.time_step_operator.time_end=0.01
                      --parser_context.rng.type=function
                      --parser_context.rng.expression=x,y:0
                      --model.time_step_operator.linear_solver.layout.writer.svg.path=mitchell_schaefer_pattern.svg)
    add_test(
      NAME      mitchell-schaefer-pattern
      COMMAND   ${CMAKE_COMMAND} -E compare_files
                mitchell_schaefer_pattern.svg
                ${CMAKE_CURRENT_SOURCE_DIR}/data/pattern/mitchell_schaefer_pattern.svg)
    set_tests_properties(mitchell-schaefer-pattern PROPERTIES DEPENDS mitchell-schaefer)

    foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
      dune_add_test(
        NAME            poisson-${dim}D
        TIMEOUT         50
        LABELS          system
        TARGET          Dune::Copasi::Executable
        CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/poisson.ini
                        --grid.dimension=${dim}
                        --parser_context.dim.value=${dim})
    endforeach()

  endif()

  foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
    dune_add_test(
      NAME            exp-${dim}D
      TIMEOUT         10
      LABELS          system
      TARGET          Dune::Copasi::Executable
      CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/exp.ini
                      --grid.dimension=${dim})

    set(grid_extensions "2")
    set(grid_origin "-1")
    foreach(i RANGE 2 ${dim})
      set(grid_extensions "${grid_extensions} 2")
      set(grid_origin "${grid_origin} -1")
    endforeach()

    foreach(parser ${DUNE_COPASI_PARSER_LIST})
      dune_add_test(
        NAME            gauss-${dim}D-${parser}
        TIMEOUT         70
        LABELS          system
        TARGET          Dune::Copasi::Executable
        CMD_ARGS        --config=${CMAKE_CURRENT_SOURCE_DIR}/gauss.ini
                        --grid.dimension=${dim}
                        --grid.extensions=${grid_extensions}
                        --grid.origin=${grid_origin}
                        --model.parser_type=${parser})
    endforeach()
  endforeach()

endif()
