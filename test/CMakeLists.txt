enable_testing()

dune_symlink_to_source_files(FILES "grids" DESTINATION "." )
dune_symlink_to_source_files(FILES "reference" DESTINATION "." )

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

if (dune-testtools_FOUND)

  add_system_test_per_target(
    TARGET          dune_copasi
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_cell.mini)

  add_system_test_per_target(
    TARGET          dune_copasi
    SCRIPT          dune_vtkcompare.py
    INIFILE         test_exp.mini)

  # make sure dune_copasi is built by default
  set_property(TARGET dune_copasi PROPERTY EXCLUDE_FROM_ALL 0)

endif()