add_executable(dune_copasi_md dune_copasi_md.cc)
target_link_libraries(dune_copasi_md PRIVATE dune_copasi_md_lib-static)

add_executable(dune_copasi_sd EXCLUDE_FROM_ALL dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd PRIVATE dune_copasi_sd_lib-static)

add_custom_target(dune_copasi ALL DEPENDS dune_copasi_md)

install(
TARGETS dune_copasi_md
RUNTIME
  COMPONENT Executable
  DESTINATION "${CMAKE_INSTALL_BINDIR}/")

install(
TARGETS dune_copasi_sd
RUNTIME
OPTIONAL
  COMPONENT Executable
  DESTINATION "${CMAKE_INSTALL_BINDIR}/")