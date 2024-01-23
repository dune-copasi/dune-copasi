install(FILES bit_flags.hh
              data_context.hh
              exceptions.hh
              factory.hh
              fmt_style.hh
              ostream_redirect.hh
              ostream_to_spdlog.hh
              stepper.hh
              tiff_file.hh
              tiff_grayscale.hh
        COMPONENT Development
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/common")
