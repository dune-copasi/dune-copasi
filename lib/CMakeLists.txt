include(GNUInstallDirs)

target_include_directories(dune-copasi
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}"
)

# SingleDomain library
add_library(singldomain-lib
  EXCLUDE_FROM_ALL
    dune_copasi_sd_fv.cc
    dune_copasi_sd_cg.cc
    dune_copasi_sd_fv_cg.cc
)
target_link_directories(singldomain-lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>
)
target_link_libraries(singldomain-lib PRIVATE dune-copasi)
set_target_properties(singldomain-lib PROPERTIES LIBRARY_OUTPUT_NAME dunecopasisd)
set_target_properties(singldomain-lib PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasisd)

if(DUNE_COPASI_SD_LIBRARY)
  set_target_properties(singldomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi INTERFACE singldomain-lib)
  target_compile_definitions(dune-copasi INTERFACE DUNE_COPASI_SD_LIBRARY)

  install(TARGETS singldomain-lib
    EXPORT dune-copasi-targets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      COMPONENT Libraries
      NAMELINK_COMPONENT Development
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      COMPONENT Development
  )
endif()

# MultiDomain library
add_library(multidomain-lib
  EXCLUDE_FROM_ALL
    dune_copasi_md_fv.cc
    dune_copasi_md_cg.cc
    dune_copasi_md_fv_cg.cc
)
target_link_directories(multidomain-lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>
)
target_link_libraries(multidomain-lib PRIVATE dune-copasi)
set_target_properties(multidomain-lib PROPERTIES LIBRARY_OUTPUT_NAME dunecopasimd)
set_target_properties(multidomain-lib PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasimd)

if(DUNE_COPASI_MD_LIBRARY)
  set_target_properties(multidomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi INTERFACE multidomain-lib)
  target_compile_definitions(dune-copasi INTERFACE DUNE_COPASI_MD_LIBRARY)

  install(TARGETS multidomain-lib
    EXPORT dune-copasi-targets
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      COMPONENT Libraries
      NAMELINK_COMPONENT Development
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      COMPONENT Development
  )
endif()
