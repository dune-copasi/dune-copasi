# dune-copasi target
add_custom_target(dune_copasi ALL)

# Single domain options
set(DUNE_COPASI_SD_EXECUTABLE_LIBRARY "static"
    CACHE STRING "Library to build singldomain executable [static/shared]"
)

set(DUNE_COPASI_SD_EXECUTABLE OFF
    CACHE BOOL "Build singledomain executable by default"
)

# Single domain target
add_executable(dune_copasi_sd EXCLUDE_FROM_ALL dune_copasi_sd.cc)
target_link_libraries(dune_copasi_sd PRIVATE dune_copasi_sd_lib-${DUNE_COPASI_SD_EXECUTABLE_LIBRARY})
if(DUNE_COPASI_SD_EXECUTABLE)
  add_dependencies(dune_copasi dune_copasi_sd)
endif()

# Multiple domain options
set(DUNE_COPASI_MD_EXECUTABLE_LIBRARY "static"
    CACHE STRING "Library to build multidomain executable [static/shared]"
)

set(DUNE_COPASI_MD_EXECUTABLE ON
    CACHE BOOL "Build multidomain executable by default"
)

# Multiple domain target
add_executable(dune_copasi_md EXCLUDE_FROM_ALL dune_copasi_md.cc)
target_link_libraries(dune_copasi_md PRIVATE dune_copasi_md_lib-${DUNE_COPASI_MD_EXECUTABLE_LIBRARY})
if(DUNE_COPASI_MD_EXECUTABLE)
  add_dependencies(dune_copasi dune_copasi_md)
endif()

# install executables if they exist
install(
  TARGETS
    dune_copasi_md
    dune_copasi_sd
  EXPORT dune_copasi-export
  COMPONENT Runtime
  OPTIONAL
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}/"
)