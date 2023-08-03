
dune_add_library(dune-copasi-parser STATIC EXPORT_NAME Copasi::Parser)

add_dune_mpi_flags(dune-copasi-parser)

find_package(SymEngine 0.8.0)
find_package(ExprTk)

if(NOT (SymEngine_FOUND OR ExprTk_FOUND))
  set(MUPARSER_RQUIREMENT "REQUIRED") # we need at least one parser!
endif()
find_package(muparser ${MUPARSER_RQUIREMENT})

find_package(parafields QUIET)
# ... library not installed, we provide it from git
if (NOT parafields_FOUND)
  include(FetchContent)
  message("-- Declaring ParaFields")
  FetchContent_Declare(
    parafields
    GIT_REPOSITORY https://github.com/parafields/parafields-core.git
    GIT_TAG        3d711989be8f0838834f8c0667ae9b6af5061cbe
    )

  if(NOT parafields_POPULATED)
    message("-- Populating ParaFields")
    # TODO check if FFTW and MPI are fo
    FetchContent_Populate(parafields)
    add_subdirectory("${parafields_SOURCE_DIR}" "${parafields_BINARY_DIR}")
  endif()
endif()

target_link_libraries(dune-copasi-parser
  PUBLIC
    Dune::Copasi::Common
  PRIVATE
    fmt::fmt
    parafields::parafields
    $<TARGET_NAME_IF_EXISTS:ExprTk::ExprTk>
    $<TARGET_NAME_IF_EXISTS:muparser::muparser>
    $<$<BOOL:${SymEngine_FOUND}>:${SYMENGINE_LIBRARIES}>
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
)
