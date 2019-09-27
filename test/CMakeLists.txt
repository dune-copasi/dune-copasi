enable_testing()

dune_symlink_to_source_files(FILES "grids" DESTINATION "." )
dune_symlink_to_source_files(FILES "reference" DESTINATION "." )
dune_symlink_to_source_files(FILES "data" DESTINATION "." )

dune_add_test(NAME            test_dynamic_local_finite_element
              SOURCES         test_dynamic_local_finite_element.cc
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

if (dune-testtools_FOUND)

  add_system_test_per_target(
    TARGET          dune_copasi
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_cell.mini)

  add_system_test_per_target(
    TARGET          dune_copasi
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_exp.mini)

    add_system_test_per_target(
      TARGET          dune_copasi
      SCRIPT          dune_vtkcompare.py
      INIFILE         test_tiff.mini)

  dune_add_system_test(
    SOURCE          test_jacobian.cc
    BASENAME        dune_copasi_test_jacobian
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_jacobian.mini)

  target_link_libraries(dune_copasi_test_jacobian PRIVATE dune_copasi_lib)
  add_dependencies(build_unit_tests dune_copasi_test_jacobian)
  set_property(TEST "dune_copasi_test_jacobian_analytical" PROPERTY LABELS unit)
  set_property(TEST "dune_copasi_test_jacobian_numerical" PROPERTY LABELS unit)

  # make sure dune_copasi is built by default
  set_property(TARGET dune_copasi PROPERTY EXCLUDE_FROM_ALL 0)

endif()

add_custom_target(build_system_tests DEPENDS dune_copasi)