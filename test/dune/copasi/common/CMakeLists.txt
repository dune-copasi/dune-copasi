dune_add_test(NAME            ostream_redirect
              SOURCES         ostream_redirect.cc
              LINK_LIBRARIES  Dune::Copasi::Common
              LABELS          unit)

dune_add_test(NAME            tiff_grayscale
              SOURCES         tiff_grayscale.cc
              LINK_LIBRARIES  Dune::Copasi::Common
              LABELS          unit)
