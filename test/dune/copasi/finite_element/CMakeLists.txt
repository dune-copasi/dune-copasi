
dune_add_test(NAME            dynamic_power_local_finite_element
              SOURCES         dynamic_power.cc
              LINK_LIBRARIES  Dune::Copasi::Common
              LABELS          unit)
