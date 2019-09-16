add_library(dune_copasi_lib STATIC model_base.cc)

target_link_libraries(dune_copasi_lib PUBLIC ${DUNE_LIBS} muparser::muparser)

#install headers
install(FILES copasi.hh DESTINATION
        "${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi")
