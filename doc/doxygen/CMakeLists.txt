cmake_minimum_required(VERSION 3.9)
project(dune-copasi-doxygen)

# check if this is run as stand-alone cmake project
if("${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
  # build doxygen by default
  set(DOXYGEN_REQUIRED REQUIRED)
  set(DOXYGEN_ALL ALL)

  # get version number from dune.module
  set(REGEX "^Version:[ ]*[^\\n]+")
  file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../../dune.module" DOXYGEN_PROJECT_NUMBER REGEX "^Version:[ ]*[^\\n]+")
  set(REGEX "^[ ]*Version:[ ]*(.+)[ ]*$")
  string(REGEX REPLACE ${REGEX} "\\1" DOXYGEN_PROJECT_NUMBER "${DOXYGEN_PROJECT_NUMBER}")
endif()

find_package(Doxygen ${DOXYGEN_REQUIRED})

if(Doxygen_FOUND)
  set(DOXYGEN_EXTRACT_ALL "Yes")
  set(DOXYGEN_PROJECT_NAME "Dune::Copasi")
  set(DOXYGEN_EXCLUDE_PATTERNS "*/test/*;*/doc/*;*/src/*;*/util/*;*/python/*;*/docker/*")
  set(DOXYGEN_EXCLUDE_SYMBOLS "Impl")
  doxygen_add_docs(dune-copasi-doc WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../../" ${DOXYGEN_ALL})
endif()
