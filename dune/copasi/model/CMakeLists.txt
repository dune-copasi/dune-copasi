install(FILES constraints.hh
              functor_factory_parser.hh
              functor_factory_parser.impl.hh
              functor_factory.hh
              interpolate.hh
              local_domain.hh
              make_initial.hh
              make_step_operator.hh
              model.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/model")

add_subdirectory(diffusion_reaction)
