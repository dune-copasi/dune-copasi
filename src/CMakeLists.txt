
include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_BINDIR})

add_subdirectory(dune)

########## Create main library ##########

# create main library target
add_library(dune-copasi INTERFACE)
add_library(Dune::Copasi ALIAS dune-copasi)
set_target_properties(dune-copasi PROPERTIES EXPORT_NAME Copasi)

target_link_libraries(dune-copasi
  INTERFACE
    Dune::Copasi::Common
    Dune::Copasi::Model
    Dune::Copasi::Parser
)

# target configuration for downstream projects
install(TARGETS dune-copasi EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)

########## Create main executable ##########

option(DUNE_COPASI_WITH_EXECUTABLE "Build dune-copasi executable" ON)

if (${DUNE_COPASI_WITH_EXECUTABLE})

  include(FetchContent)

  FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
  FetchContent_MakeAvailable(json)

  add_executable(dune-copasi-exec)
  add_executable(Dune::Copasi::Executable ALIAS dune-copasi-exec)
  target_sources(dune-copasi-exec PRIVATE dune_copasi.cc)
  target_link_libraries(dune-copasi-exec
    PRIVATE
      Dune::Copasi
      nlohmann_json::nlohmann_json
  )

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

  # target configuration for downstream projects
  install(TARGETS dune-copasi-exec EXPORT dune-copasi-targets RUNTIME COMPONENT Runtime)

  # install wasm executable if compiled with emscripten
  install(
    FILES $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/>${CMAKE_INSTALL_BINDIR}/dune-copasi.wasm
    TYPE BIN
    COMPONENT Runtime
    OPTIONAL)

  set(config_opts_json "${CMAKE_CURRENT_SOURCE_DIR}/../doc/docusaurus/docs/assets/config_opts.json")
  file(READ "${config_opts_json}" CONFIGURE_OPTIONS_JSON)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${config_opts_json}")

  string(LENGTH "${CONFIGURE_OPTIONS_JSON}" CONFIGURE_OPTIONS_JSON_LENGTH)
  set(CONFIGURE_OPTIONS_JSON_LITERAL "")
  set(CONFIGURE_OPTIONS_JSON_OFFSET 0)
  set(CONFIGURE_OPTIONS_JSON_CHUNK_SIZE 12000)
  while(CONFIGURE_OPTIONS_JSON_OFFSET LESS CONFIGURE_OPTIONS_JSON_LENGTH)
    math(EXPR CONFIGURE_OPTIONS_JSON_REMAINING "${CONFIGURE_OPTIONS_JSON_LENGTH} - ${CONFIGURE_OPTIONS_JSON_OFFSET}")
    set(CONFIGURE_OPTIONS_JSON_CURRENT_CHUNK_SIZE ${CONFIGURE_OPTIONS_JSON_CHUNK_SIZE})
    if(CONFIGURE_OPTIONS_JSON_REMAINING LESS CONFIGURE_OPTIONS_JSON_CHUNK_SIZE)
      set(CONFIGURE_OPTIONS_JSON_CURRENT_CHUNK_SIZE ${CONFIGURE_OPTIONS_JSON_REMAINING})
    endif()
    string(SUBSTRING
      "${CONFIGURE_OPTIONS_JSON}"
      ${CONFIGURE_OPTIONS_JSON_OFFSET}
      ${CONFIGURE_OPTIONS_JSON_CURRENT_CHUNK_SIZE}
      CONFIGURE_OPTIONS_JSON_CHUNK
    )
    string(APPEND CONFIGURE_OPTIONS_JSON_LITERAL "R\"DCOJSON(${CONFIGURE_OPTIONS_JSON_CHUNK})DCOJSON\"\n")
    math(EXPR CONFIGURE_OPTIONS_JSON_OFFSET "${CONFIGURE_OPTIONS_JSON_OFFSET} + ${CONFIGURE_OPTIONS_JSON_CURRENT_CHUNK_SIZE}")
  endwhile()

  configure_file("config_opts.hh.in" "${PROJECT_BINARY_DIR}/config_opts.hh")

  target_include_directories(dune-copasi-exec PRIVATE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  )
endif()
