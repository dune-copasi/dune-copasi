install(FILES
        factory.hh
        local_equations.hh
        local_operator.hh
        make_step_operator.hh
        model_multi_compartment_traits.hh
        model_multi_compartment.hh
        model_multi_compartment.impl.hh
        model_single_compartment.impl.hh
        model_single_compartment.hh
        model_single_compartment_traits.hh
        reduce.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/model/diffusion_reaction")
