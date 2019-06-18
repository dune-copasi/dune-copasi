add_executable("dune-copasi-exec" dune-copasi.cc)

dune_target_link_libraries("dune-copasi-exec" ${DUNE_LIBS})
target_link_libraries("dune-copasi-exec" dune_copasi)