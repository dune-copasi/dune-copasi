add_executable(dune_copasi_md dune_copasi_md.cc)
target_link_libraries(dune_copasi_md PRIVATE dune_copasi_lib)

add_executable(dune_copasi_sd dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd PRIVATE dune_copasi_lib)

add_custom_target(dune_copasi ALL DEPENDS dune_copasi_sd dune_copasi_md)

install(TARGETS dune_copasi_md RUNTIME
        DESTINATION "${CMAKE_INSTALL_BINDIR}/")

install(TARGETS dune_copasi_sd RUNTIME
        DESTINATION "${CMAKE_INSTALL_BINDIR}/")