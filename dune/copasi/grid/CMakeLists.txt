install(FILES boundary_entity_mapper.hh
              cell_data_parser.hh
              cell_data.hh
              has_single_geometry_type.hh
              make_multi_domain_grid.hh
              mark_stripes.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/grid")
