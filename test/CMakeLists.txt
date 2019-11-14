enable_testing()

dune_symlink_to_source_files(FILES "grids" DESTINATION "." )
dune_symlink_to_source_files(FILES "reference" DESTINATION "." )
dune_symlink_to_source_files(FILES "data" DESTINATION "." )

dune_add_test(NAME            test_dynamic_power_local_finite_element
              SOURCES         test_dynamic_power_local_finite_element.cc
              LABELS          unit)

dune_add_test(NAME            test_concepts_pdelab
              SOURCES         test_concepts_pdelab.cc
              LABELS          unit)

dune_add_test(NAME            test_concepts_grid
              SOURCES         test_concepts_grid.cc
              LABELS          unit)

dune_add_test(NAME            test_concepts_typetree
              SOURCES         test_concepts_typetree.cc
              LABELS          unit)

dune_add_test(NAME            test_tiff_grayscale
              SOURCES         test_tiff_grayscale.cc
              LINK_LIBRARIES  TIFF::TIFF
              LABELS          unit)

#dune_add_test(NAME            test_parameter_tree_check
#              SOURCES         test_parameter_tree_check.cc
#              LABELS          unit)

add_custom_target(build_system_tests DEPENDS dune_copasi)
add_executable(dune_copasi_sd_compare EXCLUDE_FROM_ALL dune_copasi_sd_compare.cc)
add_executable(dune_copasi_md_compare EXCLUDE_FROM_ALL dune_copasi_md_compare.cc)
target_link_libraries(dune_copasi_sd_compare PRIVATE dune_copasi_sd_lib)
target_link_libraries(dune_copasi_md_compare PRIVATE dune_copasi_md_lib)
add_dependencies(build_system_tests dune_copasi_sd_compare dune_copasi_md_compare)

if (dune-testtools_FOUND)

  add_system_test_per_target(
    TARGET          dune_copasi_sd_compare
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_exp.mini)

  add_system_test_per_target(
    TARGET          dune_copasi_md # todo: use md compare
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_exp.mini)

  add_system_test_per_target(
    TARGET          dune_copasi_sd_compare
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_gauss.mini)

  add_system_test_per_target(
    TARGET          dune_copasi_md # todo: use md compare
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_gauss.mini)

  add_system_test_per_target(
    TARGET          dune_copasi_sd
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_tiff.mini)

  add_system_test_per_target(
    TARGET          dune_copasi_md
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_tiff.mini)

  add_system_test_per_target(
    TARGET          dune_copasi_md
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_cell.mini)

  dune_add_system_test(
    SOURCE          test_jacobian.cc
    BASENAME        dune_copasi_test_jacobian
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_jacobian.mini)

  target_link_libraries(dune_copasi_test_jacobian PRIVATE dune_copasi_lib)
  add_dependencies(build_unit_tests dune_copasi_test_jacobian)
  set_property(TEST "dune_copasi_test_jacobian_analytical" PROPERTY LABELS unit)
  set_property(TEST "dune_copasi_test_jacobian_numerical" PROPERTY LABELS unit)

  # make sure dune_copasi executables are built by default
  set_property(TARGET dune_copasi_sd PROPERTY EXCLUDE_FROM_ALL 0)
  set_property(TARGET dune_copasi_md PROPERTY EXCLUDE_FROM_ALL 0)
endif()