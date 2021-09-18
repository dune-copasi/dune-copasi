install(FILES bit_flags.hh
              enum.hh
              filesystem.hh
              local_basis_cache.hh
              parser_to_grid_function.hh
              stepper.hh
              tiff_file.hh
              tiff_grayscale.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/common")
