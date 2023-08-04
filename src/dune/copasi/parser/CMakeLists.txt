
find_package(SymEngine 0.8.0)
if(${SymEngine_FOUND})
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}\nfind_package(SymEngine 0.8.0 REQUIRED)")
endif()

find_package(ExprTk)
if(${ExprTk_FOUND})
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}\nfind_package(ExprTk REQUIRED)")
endif()

if(NOT (SymEngine_FOUND OR ExprTk_FOUND))
  set(MUPARSER_REQUIREMENT "REQUIRED") # we need at least one parser!
endif()
find_package(muparser ${MUPARSER_REQUIREMENT})
if(${muparser_FOUND})
  set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}\nfind_package(muparser REQUIRED)")
  set(DUNE_COPASI_MUPARSER_MAX_FUNCTIONS 500 CACHE STRING "Maximum number of functions that muparser can store in its context")
endif()

option(DUNE_COPASI_USE_PARAFIELDS "Use parafields (random field generator) for the dune-copasi parsers" ON)
if(${DUNE_COPASI_USE_PARAFIELDS})
  find_package(parafields QUIET)
  # ... library not installed, we provide it from git
  if (NOT parafields_FOUND)
    include(FetchContent)
    message("-- Declaring ParaFields")
    FetchContent_Declare(
      parafields
      GIT_REPOSITORY https://github.com/parafields/parafields-core.git
      GIT_TAG        1e276ec12e23dbbf8f0858dae12076e94e5fae08
    )

    if(NOT parafields_POPULATED)
      message("-- Populating ParaFields")
      FetchContent_Populate(parafields)
      add_subdirectory("${parafields_SOURCE_DIR}" "${parafields_BINARY_DIR}")
      set(parafields_FOUND ON)
    endif()
  endif()
  if(${parafields_FOUND})
    set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}\nfind_package(parafields REQUIRED)")
  endif()
endif()

set(DUNE_CUSTOM_PKG_CONFIG_SECTION "${DUNE_CUSTOM_PKG_CONFIG_SECTION}" PARENT_SCOPE)

add_library(dune-copasi-parser STATIC)
add_library(Dune::Copasi::Parser ALIAS dune-copasi-parser)
set_target_properties(dune-copasi-parser PROPERTIES EXPORT_NAME Copasi::Parser)

add_dune_mpi_flags(dune-copasi-parser)

target_link_libraries(dune-copasi-parser
  PUBLIC
    Dune::Copasi::Common
  PRIVATE
    $<TARGET_NAME_IF_EXISTS:parafields::parafields>
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

target_compile_definitions(dune-copasi-parser
  PRIVATE
    $<$<TARGET_EXISTS:muparser::muparser>:DUNE_COPASI_MUPARSER_MAX_FUNCTIONS=${DUNE_COPASI_MUPARSER_MAX_FUNCTIONS}>
)

# target configuration for downstream projects
install(TARGETS dune-copasi-parser EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
