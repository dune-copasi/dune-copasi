cmake_minimum_required(VERSION 3.13)
project(dune-copasi-tests CXX)

if(NOT TARGET dune-copasi::dune-copasi)
  find_package(dune-copasi REQUIRED IMPORTED COMPONENTS Libraries Development Runtime)
endif()

if(NOT dune-testtools)
  find_package(dune-testtools)
endif()

if (dune-testtools_FOUND)
  list(APPEND CMAKE_MODULE_PATH "${dune-testtools_MODULE_PATH}")
  include(DuneTesttoolsMacros)
  # check that testtools python module is installed
  file(WRITE
    "${CMAKE_CURRENT_BINARY_DIR}/check_testtools.py"
    "import dune.testtools"
  )
  dune_execute_process(
    COMMAND
      ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/check_testtools.py"
    ERROR_MESSAGE
      "Python module for dune.testtools is not available."
      "Run `make install_python` on dune-testtools build directory"
      "Ensure to pass the same DUNE_PYTHON_VIRTUALENV_PATH flag as for dune-testtools"
  )

endif()

enable_testing()

dune_symlink_to_source_files(FILES "grids" DESTINATION "." )
dune_symlink_to_source_files(FILES "reference" DESTINATION "." )
dune_symlink_to_source_files(FILES "data" DESTINATION "." )

dune_add_test(NAME            test_dynamic_power_local_finite_element
              SOURCES         test_dynamic_power_local_finite_element.cc
              LINK_LIBRARIES  dune-copasi::dune-copasi
              LABELS          unit)

dune_add_test(NAME            test_concepts_pdelab
              SOURCES         test_concepts_pdelab.cc
              LINK_LIBRARIES  dune-copasi::dune-copasi
              LABELS          unit)

dune_add_test(NAME            test_concepts_grid
              SOURCES         test_concepts_grid.cc
              LINK_LIBRARIES  dune-copasi::dune-copasi
              LABELS          unit)

dune_add_test(NAME            test_concepts_typetree
              SOURCES         test_concepts_typetree.cc
              LINK_LIBRARIES  dune-copasi::dune-copasi
              LABELS          unit)

dune_add_test(NAME            test_tiff_grayscale
              SOURCES         test_tiff_grayscale.cc
              LINK_LIBRARIES  dune-copasi::dune-copasi
              LABELS          unit)

# # SingleDomain tests
add_executable(dune-copasi-sd-compare EXCLUDE_FROM_ALL dune_copasi_sd_compare.cc)

target_link_libraries(dune-copasi-sd-compare PRIVATE dune-copasi::dune-copasi)

add_custom_target(build_system_tests)
add_dependencies(build_system_tests dune-copasi-sd-compare)

if (dune-testtools_FOUND)

  add_system_test_per_target(
    TARGET          dune-copasi-sd-compare
    INIFILE         test_exp.mini)

  add_system_test_per_target(
    TARGET          dune-copasi-sd-compare
    INIFILE         test_initial.mini)

  add_system_test_per_target(
    TARGET          dune-copasi-sd-compare
    INIFILE         test_gauss.mini)

  if(dune-copasi::singledomain-exec)
    add_system_test_per_target(
      TARGET          dune-copasi::singledomain-exec
      SCRIPT          dune_vtkcompare.py
      INIFILE         test_tiff.mini)

    add_system_test_per_target(
      TARGET          dune-copasi::singledomain-exec
      INIFILE         test_brusselator.mini)

    add_system_test_per_target(
      TARGET          dune-copasi::singledomain-exec
      INIFILE         test_no_writer.mini)
  endif()
endif()

# MultiDomain tests
add_executable(dune-copasi-md-compare EXCLUDE_FROM_ALL dune_copasi_md_compare.cc)
add_executable(dune-copasi-md-jacobian EXCLUDE_FROM_ALL test_jacobian.cc)

target_link_libraries(dune-copasi-md-compare PRIVATE dune-copasi::dune-copasi)
target_link_libraries(dune-copasi-md-jacobian PRIVATE dune-copasi::dune-copasi)

add_dependencies(build_system_tests
  dune-copasi-md-compare
)

add_dependencies(build_unit_tests
  dune-copasi-md-jacobian
)

if (dune-testtools_FOUND)

  add_system_test_per_target(
    TARGET          dune-copasi-md-compare
    INIFILE         test_exp.mini)


  add_system_test_per_target(
    TARGET          dune-copasi-md-compare
    INIFILE         test_initial.mini)

  add_system_test_per_target(
    TARGET          dune-copasi-md-compare
    INIFILE         test_gauss.mini)

  add_system_test_per_target(
    TARGET          dune-copasi-md-compare
    INIFILE         test_gauss_mixed_geometry.mini)

  if(dune-copasi::multidomain-exec)
    add_system_test_per_target(
      TARGET          dune-copasi::multidomain-exec
      SCRIPT          dune_vtkcompare.py
      INIFILE         test_tiff.mini)

    add_system_test_per_target(
      TARGET          dune-copasi::multidomain-exec
      SCRIPT          dune_vtkcompare.py
      INIFILE         test_cell.mini)

    add_system_test_per_target(
      TARGET          dune-copasi::multidomain-exec
      INIFILE         test_no_writer.mini)
  endif()

  add_system_test_per_target(
    TARGET          dune-copasi-md-jacobian
    INIFILE         test_cell.mini)

  set_property(TEST "dune-copasi-md-jacobian_test_cell" PROPERTY LABELS unit)
  set_property(TEST "dune-copasi-md-jacobian_test_cell" PROPERTY LABELS unit)

  add_system_test_per_target(
    TARGET          dune-copasi-md-jacobian
    INIFILE         test_NFkappaB.mini)

  set_property(TEST "dune-copasi-md-jacobian_test_NFkappaB" PROPERTY LABELS unit)
  set_property(TEST "dune-copasi-md-jacobian_test_NFkappaB" PROPERTY LABELS unit)

endif()
