cmake_minimum_required(VERSION 3.13)
project(dune-copasi-tests CXX)

# check if this project is called as its own cmake project (i.e. not as subdirectory)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # ...if so, find dune-copasi as if this were downstream project
  find_package(dune-copasi REQUIRED IMPORTED COMPONENTS Libraries Development Runtime)
endif()

add_subdirectory(dune)

# if(NOT dune-testtools)
#   find_package(dune-testtools)
# endif()

# if (dune-testtools_FOUND)
#   list(APPEND CMAKE_MODULE_PATH "${dune-testtools_MODULE_PATH}")
#   include(DuneTesttoolsMacros)
#   # check that testtools python module is installed
#   file(WRITE
#     "${CMAKE_CURRENT_BINARY_DIR}/check_testtools.py"
#     "import dune.testtools"
#   )
#   dune_execute_process(
#     COMMAND
#       ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/check_testtools.py"
#     ERROR_MESSAGE
#       "Python module for dune.testtools is not available."
#       "Run `make install_python` on dune-testtools build directory"
#       "Ensure to pass the same DUNE_PYTHON_VIRTUALENV_PATH flag as for dune-testtools"
#   )
# endif()

enable_testing()

# dune_symlink_to_source_files(FILES "grids" DESTINATION "." )
# dune_symlink_to_source_files(FILES "reference" DESTINATION "." )
# dune_symlink_to_source_files(FILES "data" DESTINATION "." )

# dune_add_test(NAME            test_dynamic_power_local_finite_element
#               SOURCES         test_dynamic_power_local_finite_element.cc
#               LINK_LIBRARIES  Dune::Copasi::Common
#               LABELS          unit)

# dune_add_test(NAME            test_event_stepper
#               SOURCES         test_event_stepper.cc
#               LINK_LIBRARIES  Dune::Copasi::Common
#               LABELS          unit)

# # # # SingleDomain tests
# # add_executable(dune-copasi-sd-compare EXCLUDE_FROM_ALL dune_copasi_sd_compare.cc)

# # target_link_libraries(dune-copasi-sd-compare PRIVATE Dune::Copasi)

# add_custom_target(build_system_tests)
# # add_dependencies(build_system_tests dune-copasi-sd-compare)

# if (dune-testtools_FOUND)

#   # add_system_test_per_target(
#   #   TARGET          dune-copasi-sd-compare
#   #   INIFILE         test_initial.mini)

#   # add_system_test_per_target(
#   #   TARGET          dune-copasi-sd-compare
#   #   INIFILE         test_gauss.mini)

#   if(TARGET Dune::Copasi::Executable2D)

#   #   add_system_test_per_target(
#   #     TARGET          Dune::Copasi::Executable2D
#   #     SCRIPT          dune_vtkcompare.py
#   #     INIFILE         test_tiff.mini)

#     add_system_test_per_target(
#       TARGET          Dune::Copasi::Executable2D
#       INIFILE         test_brusselator.mini)

#   #   add_system_test_per_target(
#   #     TARGET          Dune::Copasi::Executable2D
#   #     INIFILE         test_time_snap.mini)

#   #   add_system_test_per_target(
#   #     TARGET          Dune::Copasi::Executable2D
#   #     INIFILE         test_no_writer.mini)
#   endif()
# endif()

# dune_add_test(
#   TARGET          Dune::Copasi::Executable2D
#   CMD_ARGS        ${CMAKE_CURRENT_SOURCE_DIR}/test_exp.ini)

# if (dune-testtools_FOUND)

#   # add_system_test_per_target(
#   #   TARGET          Dune::Copasi::Executable2D
#   #   INIFILE         test_exp.mini)


#   add_system_test_per_target(
#     TARGET          Dune::Copasi::Executable2D
#     INIFILE         test_initial.mini)

#   add_system_test_per_target(
#     TARGET          Dune::Copasi::Executable2D
#     INIFILE         test_gauss.mini)

#   if(TARGET Dune::Copasi::Executable2D)
#     # add_system_test_per_target(
#     #   TARGET          Dune::Copasi::Executable2D
#     #   SCRIPT          dune_vtkcompare.py
#     #   INIFILE         test_cell.mini)

#     add_system_test_per_target(
#       TARGET          Dune::Copasi::Executable2D
#       INIFILE         test_membrane.mini)
#   endif()

#   add_system_test_per_target(
#     TARGET          Dune::Copasi::Executable2D
#     INIFILE         test_cell.mini)

#   # set_property(TEST "dune-copasi-md-jacobian_test_cell" PROPERTY LABELS unit)
#   # set_property(TEST "dune-copasi-md-jacobian_test_cell" PROPERTY LABELS unit)

#   add_system_test_per_target(
#     TARGET          Dune::Copasi::Executable2D
#     INIFILE         test_NFkappaB.mini)

#   # set_property(TEST "dune-copasi-md-jacobian_test_NFkappaB" PROPERTY LABELS unit)
#   # set_property(TEST "dune-copasi-md-jacobian_test_NFkappaB" PROPERTY LABELS unit)

# endif()

if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  set_target_properties(build_unit_tests build_system_tests PROPERTIES EXCLUDE_FROM_ALL OFF)
endif()
