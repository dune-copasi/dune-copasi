add_subdirectory(diffusion_reaction)

install(FILES variadic.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/local_operator")
