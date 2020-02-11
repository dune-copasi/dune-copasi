install(FILES multidomain_gmsh_reader.hh
              multidomain_entity_transformation.hh
              has_single_geometry_type.hh
              mark_stripes.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/grid")
