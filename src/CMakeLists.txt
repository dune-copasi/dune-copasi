add_executable(dune_copasi dune_copasi.cc)
target_link_libraries(dune_copasi PRIVATE dune_copasi_lib)

install(TARGETS dune_copasi RUNTIME
        DESTINATION "${CMAKE_INSTALL_BINDIR}/")