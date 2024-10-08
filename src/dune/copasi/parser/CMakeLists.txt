
set(parser_list)
find_package(SymEngine 0.8.0)
if(${SymEngine_FOUND})
  list(APPEND parser_list SymEngine SymEngineSBML)
  add_library(Dune::Wrapper::SymEngine INTERFACE IMPORTED GLOBAL)
  target_link_libraries(Dune::Wrapper::SymEngine INTERFACE ${SYMENGINE_LIBRARIES})
  target_include_directories(Dune::Wrapper::SymEngine INTERFACE ${SYMENGINE_INCLUDE_DIRS})
endif()

find_package(ExprTk)
option(DUNE_COPASI_DISABLE_FETCH_PACKAGE_ExprTk "Disable (internet) fetch package for ExprTk" OFF)
if ((NOT ExprTk_FOUND) AND (NOT CMAKE_DISABLE_FIND_PACKAGE_ExprTk) AND (NOT DUNE_COPASI_DISABLE_FETCH_PACKAGE_ExprTk))
  include(FetchContent)
  message("-- Declaring ExprTk")
  FetchContent_Declare(
    exprtk
    GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git
    GIT_TAG        0.0.2
  )
  if(NOT exprtk_POPULATED)
    message("-- Populating ExprTk")
    FetchContent_Populate(exprtk)
    set(ExprTk_INCLUDE_DIRS "${exprtk_SOURCE_DIR}")
    set(ExprTk_FOUND ON)
  endif()
endif()
if(${ExprTk_FOUND})
  list(APPEND parser_list ExprTk)
endif()

find_package(muparser)
if(${muparser_FOUND})
  set(DUNE_COPASI_MUPARSER_MAX_FUNCTIONS 5000 CACHE STRING "Maximum number of functions that muparser can store in its context")
  list(APPEND parser_list MuParser)
endif()

find_package(parafields)
option(DUNE_COPASI_DISABLE_FETCH_PACKAGE_parafields "Disable (internet) fetch package for parafields" ON)
if((NOT parafields_FOUND) AND (NOT CMAKE_DISABLE_FIND_PACKAGE_parafields) AND (NOT DUNE_COPASI_DISABLE_FETCH_PACKAGE_parafields))
  include(FetchContent)
  message("-- Declaring ParaFields")
  FetchContent_Declare(
    parafields
    GIT_REPOSITORY https://github.com/parafields/parafields-core.git
    GIT_TAG        61b2893f4940e3bb5881bf32faee14524efb4184
  )

  if(NOT parafields_POPULATED)
    message("-- Populating ParaFields")
    FetchContent_Populate(parafields)
    set(BUILD_TESTING OFF)
    add_subdirectory("${parafields_SOURCE_DIR}" "${parafields_BINARY_DIR}")
    set(parafields_FOUND ON)
  endif()
endif()

set(DUNE_COPASI_PARSER_LIST ${parser_list} CACHE INTERNAL "List of parsers found found and used in dune-copasi-parser")

add_library(dune-copasi-parser STATIC)
add_library(Dune::Copasi::Parser ALIAS dune-copasi-parser)
set_target_properties(dune-copasi-parser PROPERTIES EXPORT_NAME Copasi::Parser)

add_dune_mpi_flags(dune-copasi-parser)

target_link_libraries(dune-copasi-parser
  PUBLIC
    Dune::Copasi::Common
  PRIVATE
    $<TARGET_NAME_IF_EXISTS:parafields::parafields>
    $<TARGET_NAME_IF_EXISTS:muparser::muparser>
    $<TARGET_NAME_IF_EXISTS:Dune::Wrapper::SymEngine>
  )

target_sources(dune-copasi-parser
  PRIVATE
    factory.cc
    context.cc
    parser.cc
    $<$<BOOL:${ExprTk_FOUND}>:${CMAKE_CURRENT_SOURCE_DIR}/exprtk.cc>
    $<$<TARGET_EXISTS:muparser::muparser>:${CMAKE_CURRENT_SOURCE_DIR}/mu.cc>
    $<$<TARGET_EXISTS:Dune::Wrapper::SymEngine>:${CMAKE_CURRENT_SOURCE_DIR}/symengine.cc>
)

target_include_directories(dune-copasi-parser
  PRIVATE
    $<$<BOOL:${ExprTk_FOUND}>:${ExprTk_INCLUDE_DIRS}>
)

target_compile_definitions(dune-copasi-parser
  PRIVATE
    $<$<TARGET_EXISTS:muparser::muparser>:DUNE_COPASI_MUPARSER_MAX_FUNCTIONS=${DUNE_COPASI_MUPARSER_MAX_FUNCTIONS}>
    $<$<TARGET_EXISTS:muparser::muparser>:DUNE_COPASI_HAVE_MUPARSER>
)

# target configuration for downstream projects
install(TARGETS dune-copasi-parser EXPORT dune-copasi-targets ARCHIVE COMPONENT Development)
