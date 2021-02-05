include(GNUInstallDirs)

target_include_directories(dune-copasi
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}"
)

unset(EXPORTED_DEVELOPMENT_TARGETS)

# SingleDomain library
add_library(singledomain-lib
  EXCLUDE_FROM_ALL
    dune_copasi_sd_fv.cc
    dune_copasi_sd_cg.cc
    dune_copasi_sd_fv_cg.cc
)
target_link_directories(singledomain-lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>
)
target_link_libraries(singledomain-lib PRIVATE dune-copasi)
set_target_properties(singledomain-lib PROPERTIES LIBRARY_OUTPUT_NAME dunecopasisd)
set_target_properties(singledomain-lib PROPERTIES ARCHIVE_OUTPUT_NAME dunecopasisd)

if(DUNE_COPASI_SD_LIBRARY)
  set_target_properties(singledomain-lib PROPERTIES EXCLUDE_FROM_ALL OFF)
  target_link_libraries(dune-copasi INTERFACE singledomain-lib)
  target_compile_definitions(dune-copasi INTERFACE DUNE_COPASI_SD_LIBRARY)
  list(APPEND EXPORTED_DEVELOPMENT_TARGETS singledomain-lib)
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
  list(APPEND EXPORTED_DEVELOPMENT_TARGETS multidomain-lib)
endif()

# Set up filesystem
add_library(filesystem INTERFACE)
list(APPEND EXPORTED_DEVELOPMENT_TARGETS filesystem)

if(Filesystem_FOUND)
  # file system is supported, link and we are done
  target_link_libraries(filesystem INTERFACE std::filesystem)
else()
  # ... else, we provide a fallback
  include(FetchContent)
  message("-- Declaring GHC Filesystem")
  FetchContent_Declare(
    ghc_filesystem
    GIT_REPOSITORY https://github.com/gulrak/filesystem
    GIT_TAG        v1.4.0
  )

  if(NOT ghc_filesystem_POPULATED)
    message("-- Populating GHC Filesystem")
    
    FetchContent_Populate(ghc_filesystem)
    add_subdirectory("${ghc_filesystem_SOURCE_DIR}" "${ghc_filesystem_BINARY_DIR}")
   
    install(DIRECTORY "${ghc_filesystem_SOURCE_DIR}/include"
      COMPONENT Development
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
  endif()

  target_include_directories(filesystem INTERFACE
    $<BUILD_INTERFACE:${ghc_filesystem_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

  target_compile_definitions(filesystem INTERFACE DUNE_COPASI_USE_FALLBACK_FILESYSTEM)
endif()

# main library consumes filesystem definitions
target_link_libraries(dune-copasi PUBLIC filesystem)

# install development targets
install(TARGETS ${EXPORTED_DEVELOPMENT_TARGETS}
  EXPORT dune-copasi-targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT Libraries
    NAMELINK_COMPONENT Development
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    COMPONENT Development
)

include(SetupDuneCopasiCPack)