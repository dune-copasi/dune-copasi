include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Global flag to cause add_library() to create shared libraries if ON" OFF)

target_include_directories(dune-copasi
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

# SingleDomain library
add_library(dune-copasi-sd-lib
  EXCLUDE_FROM_ALL
    "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_sd_fv.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_sd_cg.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_sd_fv_cg.cc"
)
if(DUNE_COPASI_SD_LIBRARY)
  set_target_properties(dune-copasi-sd-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi PUBLIC dune-copasi-sd-lib)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_SD_LIBRARY)
  install(TARGETS dune-copasi-sd-lib
    EXPORT dune-copasi-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    OPTIONAL
  )
endif()


# MultiDomain library
add_library(dune-copasi-md-lib
  EXCLUDE_FROM_ALL
    "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_md_fv.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_md_cg.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/dune_copasi_md_fv_cg.cc"
)
if(DUNE_COPASI_MD_LIBRARY)
  set_target_properties(dune-copasi-md-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi PUBLIC dune-copasi-md-lib)
  target_compile_definitions(dune-copasi PUBLIC DUNE_COPASI_MD_LIBRARY)
  install(TARGETS dune-copasi-md-lib
    EXPORT dune-copasi-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    OPTIONAL
  )
endif()
