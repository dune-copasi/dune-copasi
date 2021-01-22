include(GNUInstallDirs)

# SingleDomain executable
if(DUNE_COPASI_SD_EXECUTABLE)
  add_executable(dune_copasi_sd dune_copasi_sd.cc)
  target_link_libraries(dune_copasi_sd PRIVATE dune_copasi_sd_lib)
  add_dependencies(dune_copasi dune_copasi_sd)
  install(TARGETS dune_copasi_sd
    COMPONENT SigleDomain
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# MultiDomain executable
if(DUNE_COPASI_MD_EXECUTABLE)
  add_executable(dune_copasi_md dune_copasi_md.cc)
  target_link_libraries(dune_copasi_md PRIVATE dune_copasi_md_lib)
  add_dependencies(dune_copasi dune_copasi_md)
  install(TARGETS dune_copasi_md
    COMPONENT MultiDomain
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
