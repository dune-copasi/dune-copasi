install(FILES block_jacobi.hh
              concepts.hh
              dense_inverse.hh
              umfpack.hh
              util.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/solver/istl")

add_subdirectory(factory)
