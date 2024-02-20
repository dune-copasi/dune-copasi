dune_add_test(NAME            ostream_redirect
              SOURCES         ostream_redirect.cc
              LINK_LIBRARIES  Dune::Copasi::Common GTest::gtest_main
              LABELS          unit
              COMPILE_ONLY)
gtest_discover_tests(ostream_redirect PROPERTIES LABELS unit)

dune_add_test(NAME                  tiff_grayscale
              SOURCES               tiff_grayscale.cc
              LINK_LIBRARIES        Dune::Copasi::Common GTest::gtest_main
              LABELS                unit
              COMPILE_ONLY)
gtest_discover_tests(tiff_grayscale 
    PROPERTIES
        LABELS unit
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../../
)