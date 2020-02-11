add_library(dune_copasi_lib STATIC ../dune/copasi/model/base.cc)
target_link_libraries(dune_copasi_lib PUBLIC ${DUNE_LIBS} TIFF::TIFF muparser::muparser)
target_include_directories(dune_copasi_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(dune_copasi_md_lib-static EXCLUDE_FROM_ALL STATIC dune_copasi_md.cc)
target_link_libraries(dune_copasi_md_lib-static PUBLIC dune_copasi_lib)

add_library(dune_copasi_md_lib-shared EXCLUDE_FROM_ALL SHARED dune_copasi_md.cc)
target_link_libraries(dune_copasi_md_lib-shared PUBLIC dune_copasi_lib)

add_library(dune_copasi_sd_lib-static EXCLUDE_FROM_ALL STATIC dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd_lib-static PUBLIC dune_copasi_lib)

add_library(dune_copasi_sd_lib-shared EXCLUDE_FROM_ALL SHARED dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd_lib-shared PUBLIC dune_copasi_lib)

install(
  TARGETS dune_copasi_lib
          dune_copasi_md_lib-static
          dune_copasi_md_lib-shared
          dune_copasi_sd_lib-static
          dune_copasi_sd_lib-shared
  OPTIONAL
  EXPORT dune_copasi-export
  COMPONENT Library
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)