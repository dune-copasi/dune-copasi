install(FILES base.hh
              continuous_galerkin.hh
              finite_volume.hh
              multidomain.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/local_operator/diffusion_reaction")