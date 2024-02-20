
dune_add_test(NAME            dynamic_power_local_finite_element
              SOURCES         dynamic_power.cc
              LINK_LIBRARIES  Dune::Copasi::Common GTest::gtest_main
              LABELS          unit
              COMPILE_ONLY)
gtest_discover_tests(dynamic_power_local_finite_element PROPERTIES LABELS unit)