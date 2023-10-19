install(FILES diffusion_reaction.cc
              diffusion_reaction.hh
              multidomain_diffusion_reaction.cc
              multidomain_diffusion_reaction.hh
              state.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/model")
