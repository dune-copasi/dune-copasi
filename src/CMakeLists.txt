add_executable("dune_copasi_exec" dune-copasi.cc)

target_link_libraries("dune_copasi_exec" PRIVATE muparser::muparser dune_copasi)
# target_link_libraries(dune_copasi_exec PRIVATE )
# dune_target_link_libraries("dune-copasi-exec" ${DUNE_LIBS})

