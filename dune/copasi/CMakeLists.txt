add_library(dune_copasi STATIC
  model_base.cc)

target_link_libraries(dune_copasi PUBLIC ${DUNE_LIBS})

#install headers
install(FILES copasi.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi)
