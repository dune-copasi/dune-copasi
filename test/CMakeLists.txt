dune_add_test(NAME "test_run"  
              SOURCES   ${CMAKE_SOURCE_DIR}/src/dune-copasi.cc # TODO: use target instead of source
              CMD_ARGS  ${CMAKE_SOURCE_DIR}/doc/config.ini
              LABELS    "system")