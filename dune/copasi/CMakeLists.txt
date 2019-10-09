add_subdirectory(concepts)

install(FILES   coefficient_mapper.hh
                dynamic_local_basis.hh
                dynamic_local_coefficients.hh
                dynamic_local_finite_element.hh
                dynamic_local_interpolation.hh
                enum.hh
                gmsh_reader.hh
                grid_function_writer.hh
                local_operator_multidomain.hh
                local_operator.hh
                model_base.hh
                model_diffusion_reaction.cc
                model_diffusion_reaction.hh
                model_multidomain_diffusion_reaction.cc
                model_multidomain_diffusion_reaction.hh
                model_state.hh
                multidomain_entity_transformation.hh
                multidomain_local_finite_element_map.hh
                muparser_data_handler.hh
                pdelab_callable_adapter.hh
                pdelab_expression_adapter.hh
                tiff_grayscale.hh
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi")
