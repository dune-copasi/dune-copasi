install(FILES direct.hh
              inverse.hh
              iterative.hh
              preconditioner.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/solver/istl/factory")
