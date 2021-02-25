# SingleDomain objects
if(DUNE_COPASI_SD_LIBRARY)
  add_library(singledomain-lib
    OBJECT
    EXCLUDE_FROM_ALL
      dune_copasi_sd_fv.cc
      dune_copasi_sd_cg.cc
      dune_copasi_sd_fv_cg.cc
  )
  target_include_directories(singledomain-lib PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
  target_link_libraries(singledomain-lib PRIVATE dune-copasi-deps)
  target_sources(dune-copasi PRIVATE $<TARGET_OBJECTS:singledomain-lib>)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_MD_LIBRARY)
endif()

# MultiDomain objects
if(DUNE_COPASI_MD_LIBRARY)
  add_library(multidomain-lib
    OBJECT
    EXCLUDE_FROM_ALL
      dune_copasi_md_fv.cc
      dune_copasi_md_cg.cc
      dune_copasi_md_fv_cg.cc
  )
  target_include_directories(multidomain-lib PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
  target_link_libraries(multidomain-lib PRIVATE dune-copasi-deps)
  target_sources(dune-copasi PRIVATE $<TARGET_OBJECTS:multidomain-lib>)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_SD_LIBRARY)
endif()
