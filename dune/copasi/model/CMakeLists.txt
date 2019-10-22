install(FILES   base.hh
                multidomain_local_operator.hh
                local_operator.hh
                diffusion_reaction.cc
                diffusion_reaction.hh
                multidomain_diffusion_reaction.cc
                multidomain_diffusion_reaction.hh
                state.hh
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/model")