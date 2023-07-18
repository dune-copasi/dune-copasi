
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

set(target dune-copasi-exec)
add_executable(dune-copasi-exec)
add_executable(Dune::Copasi::Executable ALIAS dune-copasi-exec)
target_sources(dune-copasi-exec PRIVATE dune_copasi.cc)
target_link_libraries(dune-copasi-exec PRIVATE Dune::Copasi)

string(JOIN "," comp_def_dims ${DUNE_COPASI_GRID_DIMENSIONS})
target_compile_definitions(dune-copasi-exec PRIVATE DUNE_COPASI_GRID_DIMENSIONS=${comp_def_dims})

foreach(dim ${DUNE_COPASI_GRID_DIMENSIONS})
  string(JOIN "," comp_def_dims ${DUNE_COPASI_${dim}_FEM_ORDERS})
  target_compile_definitions(dune-copasi-exec PRIVATE DUNE_COPASI_${dim}_FEM_ORDERS=${comp_def_dims})
endforeach()

set_target_properties(dune-copasi-exec
  PROPERTIES
    RUNTIME_OUTPUT_NAME dune-copasi
    EXPORT_NAME Copasi::Executable)

# TODO: dune build system owns the export set...
#   install(TARGETS ${target}
#     EXPORT <dune-export-set>
#       COMPONENT Runtime
#       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

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
