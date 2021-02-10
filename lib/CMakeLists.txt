# SingleDomain objects
add_library(singledomain-lib
  OBJECT
  EXCLUDE_FROM_ALL
    dune_copasi_sd_fv.cc
    dune_copasi_sd_cg.cc
    dune_copasi_sd_fv_cg.cc
)
target_link_directories(singledomain-lib PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries(singledomain-lib PRIVATE dune-copasi-deps)

if(DUNE_COPASI_SD_LIBRARY)
  set_target_properties(singledomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_sources(dune-copasi PRIVATE $<TARGET_OBJECTS:singledomain-lib>)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_MD_LIBRARY)
endif()

# MultiDomain objects
add_library(multidomain-lib
  OBJECT
  EXCLUDE_FROM_ALL
    dune_copasi_md_fv.cc
    dune_copasi_md_cg.cc
    dune_copasi_md_fv_cg.cc
)
target_link_directories(multidomain-lib PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries(multidomain-lib PRIVATE dune-copasi-deps)

if(DUNE_COPASI_MD_LIBRARY)
  set_target_properties(multidomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_sources(dune-copasi PRIVATE $<TARGET_OBJECTS:multidomain-lib>)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_SD_LIBRARY)
endif()