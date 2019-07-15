dune_add_test(NAME            "test_local_finite_element"  
              SOURCES         "dynamic_local_finite_element.cc")

dune_add_test(NAME            "test_run"  
              SOURCES         ${CMAKE_SOURCE_DIR}/src/dune-copasi.cc # TODO: use target instead of source
              CMD_ARGS        ${CMAKE_BINARY_DIR}/doc/config.ini
              LINK_LIBRARIES  dune_copasi)