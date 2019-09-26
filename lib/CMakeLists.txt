add_library(dune_copasi_lib STATIC ../dune/copasi/model_base.cc
                                   model_multidomain_diffusion_reaction.cc)

target_link_libraries(dune_copasi_lib PUBLIC ${DUNE_LIBS} muparser::muparser)

install(TARGETS dune_copasi_lib 
        ARCHIVE
          DESTINATION "${CMAKE_INSTALL_LIBDIR}/")