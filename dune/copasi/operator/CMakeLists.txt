install(FILES adapter.hh
              line_search.hh
              operator.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/copasi/operator)

add_subdirectory(forward)
add_subdirectory(inverse)
