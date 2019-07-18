enable_testing()

dune_add_test(NAME            "test_local_finite_element"  
              SOURCES         "dynamic_local_finite_element.cc")

dune_add_test(NAME            "test_concepts_pdelab"  
              SOURCES         "test_concepts_pdelab.cc")

dune_add_test(NAME            "test_concepts_grid"  
              SOURCES         "test_concepts_grid.cc")

dune_add_test(NAME            "test_concepts_typetree"  
              SOURCES         "test_concepts_typetree.cc")