add_subdirectory(dynamic_power)

install(FILES   dynamic_power.hh
                p0.hh
                pk.hh
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/finite_element")