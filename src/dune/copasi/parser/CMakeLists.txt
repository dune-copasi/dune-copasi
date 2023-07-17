
dune_add_library(dune-copasi-parser STATIC EXPORT_NAME Copasi::Parser)

add_dune_mpi_flags(dune-copasi-parser)

target_compile_definitions(dune-copasi-parser
  PRIVATE
    HAVE_FFTW3_FLOAT=$<BOOL:${FFTW3_FLOAT_FOUND}>
    HAVE_FFTW3_DOUBLE=$<BOOL:${FFTW3_DOUBLE_FOUND}>
    HAVE_FFTW3_LONGDOUBLE=$<BOOL:${FFTW3_LONGDOUBLE_FOUND}>
)

target_link_libraries(dune-copasi-parser
  PUBLIC
    Dune::Copasi::Common
  PRIVATE
    fmt::fmt
    $<TARGET_NAME_IF_EXISTS:ExprTk::ExprTk>
    $<TARGET_NAME_IF_EXISTS:muparser::muparser>
    $<$<BOOL:${SymEngine_FOUND}>:${SYMENGINE_LIBRARIES}>
    $<$<BOOL:${FFTW3_FOUND}>:${FFTW3_LIBRARIES}>
  )

target_sources(dune-copasi-parser
  PRIVATE
    factory.cc
    context.cc
    parser.cc
    $<$<TARGET_EXISTS:ExprTk::ExprTk>:${CMAKE_CURRENT_SOURCE_DIR}/exprtk.cc>
    $<$<TARGET_EXISTS:muparser::muparser>:${CMAKE_CURRENT_SOURCE_DIR}/mu.cc>
    $<$<BOOL:${SymEngine_FOUND}>:${CMAKE_CURRENT_SOURCE_DIR}/symengine.cc>
)

target_include_directories(dune-copasi-parser
  PRIVATE
    $<$<BOOL:${SymEngine_FOUND}>:${SYMENGINE_INCLUDE_DIRS}>
    $<$<BOOL:${FFTW3_FOUND}>:${FFTW3_INCLUDE_DIR}>
)
