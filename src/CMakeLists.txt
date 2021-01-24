include(GNUInstallDirs)

target_sources(dune-copasi
  PRIVATE
    dune/copasi/model/base.cc
)

# SingleDomain executable
add_executable(dune-copasi-sd EXCLUDE_FROM_ALL "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_sd.cc")
target_link_libraries(dune-copasi-sd PRIVATE dune-copasi)
if(DUNE_COPASI_SD_EXECUTABLE)
  set_target_properties(dune-copasi-sd PROPERTIES EXCLUDE_FROM_ALL OFF)
  install(TARGETS dune-copasi-sd
  EXPORT dune-copasi-targets
    COMPONENT MultiDomain
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()

# MultiDomain executable

add_executable(dune-copasi-md EXCLUDE_FROM_ALL "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_md.cc")
target_link_libraries(dune-copasi-md PRIVATE dune-copasi)
if(DUNE_COPASI_MD_EXECUTABLE)
  set_target_properties(dune-copasi-md PROPERTIES EXCLUDE_FROM_ALL OFF)
  install(TARGETS dune-copasi-md
    EXPORT dune-copasi-targets
    COMPONENT MultiDomain
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
