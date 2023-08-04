include(GNUInstallDirs)

file(ARCHIVE_CREATE OUTPUT dune-copasi.1 PATHS dune-copasi.1)
install(FILES dune-copasi.1 DESTINATION ${CMAKE_INSTALL_MANDIR})
