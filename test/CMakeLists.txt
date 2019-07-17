file(COPY "data" DESTINATION ".")

dune_add_test(NAME            "test_local_finite_element"  
              SOURCES         "dynamic_local_finite_element.cc")

# dune_add_test(NAME            "test_run"  
#               SOURCES         ${CMAKE_SOURCE_DIR}/src/dune-copasi.cc # TODO: use target instead of source
#               CMD_ARGS        ${CMAKE_BINARY_DIR}/doc/config.ini
#               LINK_LIBRARIES  dune_copasi muparser::muparser TIFF::TIFF)

dune_add_test(NAME            "test_tiff_grayscale"  
              SOURCES         "tiff_grayscale.cc" # TODO: use target instead of source
              LINK_LIBRARIES  dune_copasi TIFF::TIFF) 