include(GNUInstallDirs)

target_sources(dune-copasi
  PRIVATE
    dune/copasi/model/base.cc
)

# SingleDomain executable
add_executable(singledomain-exec EXCLUDE_FROM_ALL "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_sd.cc")
target_link_libraries(singledomain-exec PRIVATE dune-copasi)
if(DUNE_COPASI_SD_EXECUTABLE)
  set_target_properties(singledomain-exec PROPERTIES EXCLUDE_FROM_ALL OFF)
  install(TARGETS singledomain-exec
  EXPORT dune-copasi-targets
    COMPONENT Runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
set_property(TARGET singledomain-exec PROPERTY RUNTIME_OUTPUT_NAME dune-copasi-sd)

# MultiDomain executable
add_executable(multidomain-exec EXCLUDE_FROM_ALL "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_md.cc")
target_link_libraries(multidomain-exec PRIVATE dune-copasi)
if(DUNE_COPASI_MD_EXECUTABLE)
  set_target_properties(multidomain-exec PROPERTIES EXCLUDE_FROM_ALL OFF)
  install(TARGETS multidomain-exec
    EXPORT dune-copasi-targets
    COMPONENT Runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
set_property(TARGET multidomain-exec PROPERTY RUNTIME_OUTPUT_NAME dune-copasi-md)
