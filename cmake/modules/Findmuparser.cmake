# Find the muparser library
#
# Usage:
#   find_package(muparser [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   muparser_FOUND               ... true if muparser is found on the system
#   muparser_LIBRARIES           ... full path to muparser library
#   muparser_INCLUDES            ... muparser include directory
#
# It defines the following targets:
#   muparser::muparser           ... muparser library to link against
#

find_path(muparser_INCLUDE_DIR muParserDef.h)
find_library(muparser_LIBRARY muparser)
mark_as_advanced(muparser_INCLUDE_DIR muparser_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  muparser
  REQUIRED_VARS muparser_LIBRARY muparser_INCLUDE_DIR
  VERSION_VAR muparser_VERSION_STRING)

if(muparser_FOUND
   AND NOT
       TARGET
       muparser::muparser)
  file(
    STRINGS
    "${muparser_INCLUDE_DIR}/muParserDef.h"
    muparser_version_str
    REGEX "^#define MUP_VERSION _T.*")
  string(
    REGEX
    REPLACE "^#define MUP_VERSION _T..(.*).."
    "\\1"
    muparser_VERSION_STRING
    "${muparser_version_str}")
  unset(muparser_version_str)

  set(muparser_test_path "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/CMakeTmp/muparser_test.cc")
  set(muparser_test_source
    "#include<muParser.h>
    int main() {
      mu::Parser p;
      p.DefineConst(\"pi\", 3.14);
      p.SetExpr(\"pi\");
      p.Eval();
    }")

  file(WRITE "${muparser_test_path}" "${muparser_test_source}")
  try_compile(COMPILE_RESULT_DYNAMIC "${CMAKE_CURRENT_BINARY_DIR}"
               SOURCES "${muparser_test_path}"
               CXX_STANDARD 17
               LINK_LIBRARIES ${muparser_LIBRARY}
               CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${muparser_INCLUDE_DIR})

  if(COMPILE_RESULT_DYNAMIC)
    add_library(
      muparser::muparser
      UNKNOWN
      IMPORTED)
    set_target_properties(
      muparser::muparser
      PROPERTIES IMPORTED_LOCATION ${muparser_LIBRARY}
                INTERFACE_INCLUDE_DIRECTORIES ${muparser_INCLUDE_DIR})
  else()
    file(WRITE "${muparser_test_path}" "${muparser_test_source}")
    try_compile(COMPILE_RESULT_STATIC "${CMAKE_CURRENT_BINARY_DIR}"
                SOURCES "${muparser_test_path}"
                COMPILE_DEFINITIONS -DMUPARSER_STATIC
                CXX_STANDARD 17
                LINK_LIBRARIES ${muparser_LIBRARY}
                CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${muparser_INCLUDE_DIR})

    if(NOT COMPILE_RESULT_STATIC)
      # message(FATAL_ERROR "A simple muparser test program could not be successfully compiled")
      set(muparser_FOUND OFF)
    else()

      add_library(
        muparser::muparser
        UNKNOWN
        IMPORTED)
      set_target_properties(
        muparser::muparser
        PROPERTIES IMPORTED_LOCATION ${muparser_LIBRARY}
                  INTERFACE_INCLUDE_DIRECTORIES ${muparser_INCLUDE_DIR}
                  INTERFACE_COMPILE_DEFINITIONS MUPARSER_STATIC)
    endif()
  endif()
endif()
