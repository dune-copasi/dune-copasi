# dune-copasi target
add_custom_target(dune_copasi ALL)

# Single domain target
set(DUNE_COPASI_SD_EXECUTABLE OFF
    CACHE BOOL "Build singledomain executable by default"
)

add_executable(dune_copasi_sd EXCLUDE_FROM_ALL dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd PRIVATE dune_copasi_sd_lib)
if(DUNE_COPASI_SD_EXECUTABLE)
  add_dependencies(dune_copasi dune_copasi_sd)
endif()

# Multiple domain target
set(DUNE_COPASI_MD_EXECUTABLE ON
    CACHE BOOL "Build multidomain executable by default"
)

add_executable(dune_copasi_md EXCLUDE_FROM_ALL dune_copasi_md.cc)
target_link_libraries(dune_copasi_md PRIVATE dune_copasi_md_lib)
if(DUNE_COPASI_MD_EXECUTABLE)
  add_dependencies(dune_copasi dune_copasi_md)
endif()