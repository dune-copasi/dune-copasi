add_subdirectory(dynamic_power)

install(FILES dynamic_power.hh
              local_basis_cache.hh
              p0.hh
              pk.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/finite_element")
