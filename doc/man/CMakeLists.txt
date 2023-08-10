include(GNUInstallDirs)

find_program(GZIP gzip)
string(TIMESTAMP CURRENT_DATE "%d %B %Y")
configure_file(dune-copasi.1.in dune-copasi.1)
if(GZIP)
  execute_process(COMMAND ${GZIP} -f dune-copasi.1 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dune-copasi.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endif()
