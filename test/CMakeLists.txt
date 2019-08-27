enable_testing()

dune_add_test(NAME            "test_local_finite_element"  
              SOURCES         "dynamic_local_finite_element.cc"
              LABELS          "unit")

dune_add_test(NAME            "test_concepts_pdelab"  
              SOURCES         "test_concepts_pdelab.cc"
              LABELS          "unit")

dune_add_test(NAME            "test_concepts_grid"  
              SOURCES         "test_concepts_grid.cc"
              LABELS          "unit")

dune_add_test(NAME            "test_concepts_typetree"  
              SOURCES         "test_concepts_typetree.cc"
              LABELS          "unit")

dune_add_test(NAME            "test_run"  
              SOURCES         ${CMAKE_SOURCE_DIR}/src/dune-copasi.cc # TODO: use target instead of source
              CMD_ARGS        ${CMAKE_BINARY_DIR}/doc/config.ini
              LINK_LIBRARIES  dune_copasi muparser::muparser
              LABELS          "system")