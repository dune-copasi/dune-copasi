# add source files to main library
target_sources(dune-copasi PRIVATE
  dune/copasi/common/tiff_file.cc
)

# set up single-domain executable
target_sources(singledomain-exec PRIVATE dune_copasi_sd.cc)
target_link_libraries(singledomain-exec PRIVATE dune-copasi)
target_compile_definitions(singledomain-exec PUBLIC $<$<BOOL:${DUNE_COPASI_COMPILE_3D}>:DUNE_COPASI_COMPILE_3D>)
set_target_properties(singledomain-exec
  PROPERTIES
    RUNTIME_OUTPUT_NAME dune-copasi-sd
    EXPORT_NAME Copasi::SingleDomainExec
    EXCLUDE_FROM_ALL $<NOT:$<BOOL:${DUNE_COPASI_SD_EXECUTABLE}>>)
if(DUNE_COPASI_SD_EXECUTABLE)
  install(TARGETS singledomain-exec
    EXPORT dune-copasi-targets
      COMPONENT Runtime
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# set up multi-domain executable
target_sources(multidomain-exec PRIVATE dune_copasi_md.cc)
target_link_libraries(multidomain-exec PRIVATE dune-copasi)
target_compile_definitions(multidomain-exec
  PUBLIC
    $<$<BOOL:${DUNE_COPASI_COMPILE_3D}>:DUNE_COPASI_COMPILE_3D>)
set_target_properties(multidomain-exec
  PROPERTIES
  RUNTIME_OUTPUT_NAME dune-copasi-md
    EXPORT_NAME Copasi::MultiDomainExec
    EXCLUDE_FROM_ALL $<NOT:$<BOOL:${DUNE_COPASI_MD_EXECUTABLE}>>)
if(DUNE_COPASI_MD_EXECUTABLE)
  install(TARGETS multidomain-exec
    EXPORT dune-copasi-targets
      COMPONENT Runtime
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# set up membrane-multi-domain executable
target_sources(multidomain-membrane-exec PRIVATE dune_copasi_md.cc)
target_link_libraries(multidomain-membrane-exec PRIVATE dune-copasi)
target_compile_definitions(multidomain-membrane-exec
  PUBLIC
    $<$<BOOL:${DUNE_COPASI_COMPILE_3D}>:DUNE_COPASI_COMPILE_3D>
    DUNE_COPASI_HAVE_MEMBRANE_SPACE=1)
set_target_properties(multidomain-membrane-exec
  PROPERTIES
  RUNTIME_OUTPUT_NAME dune-copasi-mdm
    EXPORT_NAME Copasi::MultiDomainExec
    EXCLUDE_FROM_ALL $<NOT:$<BOOL:${DUNE_COPASI_MD_EXECUTABLE}>>)
if(DUNE_COPASI_MD_EXECUTABLE)
  install(TARGETS multidomain-membrane-exec
    EXPORT dune-copasi-targets
      COMPONENT Runtime
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

