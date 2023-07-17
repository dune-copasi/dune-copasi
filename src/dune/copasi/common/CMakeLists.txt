
dune_add_library(dune-copasi-common STATIC EXPORT_NAME Copasi::Common)

target_sources(dune-copasi-common
  PRIVATE
    tiff_file.cc
    tiff_grayscale.cc
    ostream_redirect.cc
)

target_include_directories(dune-copasi-common
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${PROJECT_BINARY_INCLUDEDIR}>
)

target_compile_definitions(dune-copasi-common
  PUBLIC
    $<$<BOOL:${DUNE_USE_FALLBACK_FILESYSTEM}>:HAVE_GHC_FILESYSTEM>
)

target_link_libraries(dune-copasi-common
  PRIVATE
    fmt::fmt
    TIFF::TIFF
  INTERFACE
    spdlog::spdlog
    ${DUNE_LIBS}
)

target_compile_features(dune-copasi-common PUBLIC cxx_std_23)

target_precompile_headers(dune-copasi-common
  PUBLIC
    $<BUILD_INTERFACE:dune/common/fvector.hh>
    $<BUILD_INTERFACE:dune/common/fmatrix.hh>
    $<BUILD_INTERFACE:dune/common/hybridutilities.hh>
    $<BUILD_INTERFACE:dune/typetree/treepath.hh>
    $<BUILD_INTERFACE:dune/grid/yaspgrid.hh>
    $<BUILD_INTERFACE:dune/pdelab/common/tree_traversal.hh>
    $<BUILD_INTERFACE:dune/pdelab/common/for_each.hh>
    $<BUILD_INTERFACE:dune/copasi/common/stepper.hh>
    $<BUILD_INTERFACE:memory>
)
