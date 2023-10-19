list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

install(FILES
    FindExprTk.cmake
    Findmuparser.cmake
  COMPONENT Development
  DESTINATION ${DUNE_INSTALL_MODULEDIR})
