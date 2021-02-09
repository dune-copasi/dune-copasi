# add dune dependencies to main target
dune_target_enable_all_packages(dune-copasi)

# add dependencies
target_link_libraries(dune-copasi PUBLIC ${DUNE_LIBS} TIFF::TIFF PkgConfig::muparser)
if(NOT ${BUILD_SHARED_LIBS})
  target_compile_definitions(dune-copasi PUBLIC MUPARSER_STATIC)
endif()
if(DUNE_COPASI_COMPILE_3D)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_COMPILE_3D)
endif()

# Set up filesystem. Where do we get it from?
if(DUNE_USE_FALLBACK_FILESYSTEM)
  # and add ghc_filesystem to our filesystem target
  target_link_libraries(dune-copasi PUBLIC ghcFilesystem::ghc_filesystem)
  target_compile_definitions(dune-copasi PUBLIC DUNE_USE_FALLBACK_FILESYSTEM)
else()
  # ... else we don't need a back up. Link agaist standard library
  target_link_libraries(dune-copasi PUBLIC std::filesystem)
endif()

# define, link and install auxiliary develpment targets for dune-copasi
unset(AUX_DEVELOPMENT_TARGETS)

# SingleDomain library
add_library(singledomain-lib
  EXCLUDE_FROM_ALL
    dune_copasi_sd_fv.cc
    dune_copasi_sd_cg.cc
    dune_copasi_sd_fv_cg.cc
)
target_link_directories(singledomain-lib PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
target_link_libraries(singledomain-lib PRIVATE dune-copasi)
set_target_properties(singledomain-lib PROPERTIES LIBRARY_OUTPUT_NAME dunecopasisd)
set_target_properties(singledomain-lib PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasisd)

if(DUNE_COPASI_SD_LIBRARY)
  set_target_properties(singledomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi INTERFACE singledomain-lib)
  target_compile_definitions(dune-copasi INTERFACE DUNE_COPASI_SD_LIBRARY)
  list(APPEND AUX_DEVELOPMENT_TARGETS singledomain-lib)
endif()

# MultiDomain library
add_library(multidomain-lib
  EXCLUDE_FROM_ALL
    dune_copasi_md_fv.cc
    dune_copasi_md_cg.cc
    dune_copasi_md_fv_cg.cc
)
target_link_directories(multidomain-lib PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
target_link_libraries(multidomain-lib PRIVATE dune-copasi)
set_target_properties(multidomain-lib PROPERTIES LIBRARY_OUTPUT_NAME dunecopasimd)
set_target_properties(multidomain-lib PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasimd)

if(DUNE_COPASI_MD_LIBRARY)
  set_target_properties(multidomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi INTERFACE multidomain-lib)
  target_compile_definitions(dune-copasi INTERFACE DUNE_COPASI_MD_LIBRARY)
  list(APPEND AUX_DEVELOPMENT_TARGETS multidomain-lib)
endif()

# install development targets
install(TARGETS ${AUX_DEVELOPMENT_TARGETS}
  EXPORT dune-copasi-targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT Libraries
    NAMELINK_COMPONENT Development
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT Development
)

include(SetupDuneCopasiCPack)
