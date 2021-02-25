# add source files to main library
target_sources(dune-copasi PRIVATE dune/copasi/model/base.cc)

# set up single-domain executable
if(DUNE_COPASI_SD_EXECUTABLE)
  add_executable(singledomain-exec dune_copasi_sd.cc)
  target_link_libraries(singledomain-exec PRIVATE dune-copasi)
  set_property(TARGET singledomain-exec PROPERTY RUNTIME_OUTPUT_NAME dune-copasi-sd)
  if(DUNE_COPASI_COMPILE_3D)
    target_compile_definitions(singledomain-exec PUBLIC DUNE_COPASI_COMPILE_3D)
  endif()
  install(TARGETS singledomain-exec
  EXPORT dune-copasi-targets
    COMPONENT Runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()

# set up multi-domain executable
if(DUNE_COPASI_MD_EXECUTABLE)
  add_executable(multidomain-exec dune_copasi_md.cc)
  target_link_libraries(multidomain-exec PRIVATE dune-copasi)
  set_property(TARGET multidomain-exec PROPERTY RUNTIME_OUTPUT_NAME dune-copasi-md)
  if(DUNE_COPASI_COMPILE_3D)
    target_compile_definitions(multidomain-exec PUBLIC DUNE_COPASI_COMPILE_3D)
  endif()
  install(TARGETS multidomain-exec
    EXPORT dune-copasi-targets
    COMPONENT Runtime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()

