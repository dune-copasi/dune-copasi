
include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_BINDIR})

add_subdirectory(dune)

target_link_libraries(dune-copasi
  INTERFACE
    Dune::Copasi::Common
    Dune::Copasi::Parser
)

foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
  set(target dune-copasi-exec-${dim}D)
  add_executable(${target})
  add_executable(Dune::Copasi::Executable${dim}D ALIAS ${target})
  target_sources(${target} PRIVATE dune_copasi.cc)
  target_link_libraries(${target} PRIVATE Dune::Copasi)
  target_compile_definitions(${target}
    PRIVATE
      DUNE_COPASI_GRID_DIMENSION=${dim}
      DUNE_COPASI_MAX_FEM_ORDER=${DUNE_COPASI_${dim}D_MAX_FEM_ORDER})
  set_target_properties(${target}
    PROPERTIES
      RUNTIME_OUTPUT_NAME dune-copasi-${dim}D
      EXPORT_NAME Copasi::Executable${dim}D)
#   install(TARGETS ${target}
#     EXPORT dune-copasi-targets
#       COMPONENT Runtime
#       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()

# generate configuration file options out of file descriptors in the "options/" folder
file(GLOB opt_files "${CMAKE_CURRENT_SOURCE_DIR}/options/*")
set(opt_header "static const std::vector<std::array<std::string, 4>> config_file_opts = \{\n")
foreach(opt_file ${opt_files})
   file(STRINGS ${opt_file} opt_content)
   # the file name is the key name
   get_filename_component(opt ${opt_file} NAME)
   # first line is the type of the option
   list(POP_FRONT opt_content opt_type)
   # second line is short description
   list(POP_FRONT opt_content opt_short_doc)
   # the rest is its description
  string(JOIN "\\n" opt_content ${opt_content})
  # append to header
  set(opt_header "${opt_header}  \{\"${opt}\", \"${opt_type}\", \"${opt_short_doc}\", \"${opt_content}\" \},\n")
endforeach()

set(opt_header "${opt_header} \}\;\n")
file(WRITE "${PROJECT_BINARY_DIR}/dune-copasi-config-file-options.hh" ${opt_header})
