add_library(dune_copasi_lib STATIC ../dune/copasi/model/base.cc)
target_link_libraries(dune_copasi_lib PUBLIC ${DUNE_LIBS} TIFF::TIFF muparser::muparser)
target_include_directories(dune_copasi_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_library(dune_copasi_md_lib EXCLUDE_FROM_ALL dune_copasi_md.cc)
target_link_libraries(dune_copasi_md_lib PUBLIC dune_copasi_lib)

add_library(dune_copasi_sd_lib EXCLUDE_FROM_ALL dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd_lib PUBLIC dune_copasi_lib)