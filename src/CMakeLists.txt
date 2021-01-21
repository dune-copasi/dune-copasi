# SingleDomain executable

add_executable(dune_copasi_sd EXCLUDE_FROM_ALL dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd PRIVATE dune_copasi_sd_lib)
if(DUNE_COPASI_SD_EXECUTABLE)
  add_dependencies(dune_copasi dune_copasi_sd)
  install(TARGETS dune_copasi_sd COMPONENT SigleDomain)
endif()

# MultiDomain executable

add_executable(dune_copasi_md EXCLUDE_FROM_ALL dune_copasi_md.cc)
target_link_libraries(dune_copasi_md PRIVATE dune_copasi_md_lib)
if(DUNE_COPASI_MD_EXECUTABLE)
  add_dependencies(dune_copasi dune_copasi_md)
  install(TARGETS dune_copasi_md COMPONENT MultiDomain)
endif()
