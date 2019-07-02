dune_add_library(dune_copasi
  model_base.cc)

dune_target_link_libraries(dune_copasi ${DUNE_LIBS})
target_link_libraries(dune_copasi muparser::muparser)

#install headers
install(FILES copasi.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi)
