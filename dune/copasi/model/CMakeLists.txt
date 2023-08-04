install(FILES diffusion_reaction_mc_traits.hh
              diffusion_reaction_mc.hh
              diffusion_reaction_mc.impl.hh
              diffusion_reaction_sc_traits.hh
              diffusion_reaction_sc.hh
              diffusion_reaction_sc.impl.hh
              factory.hh
              model.hh
              interpolate.hh
              make_initial.hh
              make_step_operator.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/model")

add_subdirectory(local_equations)
