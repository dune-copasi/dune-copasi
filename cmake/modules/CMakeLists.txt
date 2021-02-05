
install(FILES
    DuneCopasiMacros.cmake
    FindFilesystem.cmake
    Findmuparser.cmake
    SetupDuneCopasiConfig.cmake
  DESTINATION ${DUNE_INSTALL_MODULEDIR})

# find dependencies
find_package(muparser REQUIRED)
find_package(TIFF REQUIRED)
find_package(Filesystem)
