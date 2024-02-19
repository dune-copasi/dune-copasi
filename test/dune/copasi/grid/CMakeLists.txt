dune_add_test(NAME            cell_data
              SOURCES         cell_data.cc
              LINK_LIBRARIES  Dune::Copasi::Common GTest::gtest
              LABELS          unit
              COMPILE_ONLY)
gtest_discover_tests(cell_data)