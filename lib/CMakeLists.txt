option(BUILD_SHARED_LIBS "Global flag to cause add_library() to create shared libraries if on" OFF)
option(DUNE_COPASI_COMPILE_3D "Compile 3D cases in libraries and executables" OFF)

add_library(dune_copasi_lib STATIC ../dune/copasi/model/base.cc)

target_link_libraries(dune_copasi_lib PUBLIC ${DUNE_LIBS} TIFF::TIFF muparser::muparser)
target_include_directories(dune_copasi_lib INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if (DUNE_COPASI_COMPILE_3D)
  target_compile_definitions(dune_copasi_lib PUBLIC DUNE_COPASI_COMPILE_3D)
endif()

if(USE_FALLBACK_FILESYSTEM)
  target_compile_definitions(dune_copasi_lib PUBLIC DUNE_COPASI_USE_FALLBACK_FILESYSTEM)
  target_link_libraries(dune_copasi_lib PUBLIC ghc_filesystem)
elseif(Filesystem_FOUND)
  target_link_libraries(dune_copasi_lib PUBLIC std::filesystem)
endif()

add_library(dune_copasi_md_lib_fv EXCLUDE_FROM_ALL dune_copasi_md_fv.cc)
target_link_libraries(dune_copasi_md_lib_fv PUBLIC dune_copasi_lib)

add_library(dune_copasi_md_lib_cg EXCLUDE_FROM_ALL dune_copasi_md_cg.cc)
target_link_libraries(dune_copasi_md_lib_cg PUBLIC dune_copasi_lib)

add_library(dune_copasi_md_lib_fv_cg EXCLUDE_FROM_ALL dune_copasi_md_fv_cg.cc)
target_link_libraries(dune_copasi_md_lib_fv_cg PUBLIC dune_copasi_lib)

add_library(dune_copasi_md_lib INTERFACE)
target_link_libraries(dune_copasi_md_lib INTERFACE
  dune_copasi_md_lib_fv
  dune_copasi_md_lib_cg
  dune_copasi_md_lib_fv_cg
)

add_library(dune_copasi_sd_lib_fv EXCLUDE_FROM_ALL dune_copasi_sd_fv.cc)
target_link_libraries(dune_copasi_sd_lib_fv PUBLIC dune_copasi_lib)

add_library(dune_copasi_sd_lib_cg EXCLUDE_FROM_ALL dune_copasi_sd_cg.cc)
target_link_libraries(dune_copasi_sd_lib_cg PUBLIC dune_copasi_lib)

add_library(dune_copasi_sd_lib_fv_cg EXCLUDE_FROM_ALL dune_copasi_sd_fv_cg.cc)
target_link_libraries(dune_copasi_sd_lib_fv_cg PUBLIC dune_copasi_lib)

add_library(dune_copasi_sd_lib INTERFACE)
target_link_libraries(dune_copasi_sd_lib INTERFACE
  dune_copasi_sd_lib_fv
  dune_copasi_sd_lib_cg
  dune_copasi_sd_lib_fv_cg
)

install(TARGETS dune_copasi_lib
  dune_copasi_md_lib_fv
  dune_copasi_md_lib_cg
  dune_copasi_md_lib_fv_cg
  dune_copasi_sd_lib_fv
  dune_copasi_sd_lib_cg
  dune_copasi_sd_lib_fv_cg
  COMPONENT Development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
