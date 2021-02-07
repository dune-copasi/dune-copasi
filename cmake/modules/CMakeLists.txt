list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

install(FILES
    FindFilesystem.cmake
    Findmuparser.cmake
  DESTINATION ${DUNE_INSTALL_MODULEDIR})