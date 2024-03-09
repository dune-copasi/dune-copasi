
# find specific dune-copasi dependencies
find_package(TIFF REQUIRED)
find_package(spdlog REQUIRED)
find_package(fmt REQUIRED)

# The function2 library provides an alternative to std::move_only_function
# (which otherwise requires c++23 in libstdc++ and is not yet implemented in libc++)
include(FetchContent)
message("-- Declaring function2")
FetchContent_Declare(
        function2
        GIT_REPOSITORY https://github.com/Naios/function2
        GIT_TAG        4.2.3
)
if(NOT function2_POPULATED)
  message("-- Populating function2")
  FetchContent_Populate(function2)
  add_library(function2 INTERFACE)
  install(DIRECTORY ${function2_SOURCE_DIR}/include/function2 DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/external)
  target_include_directories(function2
    INTERFACE
      $<BUILD_INTERFACE:${function2_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/dune/external>)
  install(TARGETS function2 EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
endif()

add_library(dune-copasi-common STATIC)
add_library(Dune::Copasi::Common ALIAS dune-copasi-common)
set_target_properties(dune-copasi-common PROPERTIES EXPORT_NAME Copasi::Common)

option(DUNE_COPASI_DISABLE_FMT_STYLE "Build dune copasi without fmt style" OFF)

target_sources(dune-copasi-common
  PRIVATE
    axis_names.cc
    tiff_file.cc
    tiff_grayscale.cc
    ostream_redirect.cc
)

target_include_directories(dune-copasi-common
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(dune-copasi-common
  PUBLIC
    function2
  PRIVATE
    TIFF::TIFF
  INTERFACE
    fmt::fmt
    spdlog::spdlog
    ${DUNE_LIBS}
)

foreach(module ${ALL_DEPENDENCIES})
  target_include_directories(dune-copasi-common PUBLIC ${${module}_INCLUDE_DIRS})
endforeach()

target_compile_definitions(dune-copasi-common PUBLIC DUNE_COPASI_DISABLE_FMT_STYLE=${DUNE_COPASI_DISABLE_FMT_STYLE})
target_compile_features(dune-copasi-common PUBLIC cxx_std_20)

# target configuration for downstream projects
install(TARGETS dune-copasi-common EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
