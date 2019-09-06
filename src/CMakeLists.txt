add_executable(dune_copasi dune-copasi.cc)
target_link_libraries(dune_copasi PRIVATE muparser::muparser dune_copasi_lib)

# target_link_libraries(dune_copasi_exec PRIVATE )
# dune_target_link_libraries("dune-copasi-exec" ${DUNE_LIBS})

