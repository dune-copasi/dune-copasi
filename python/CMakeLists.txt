configure_file("dune/copasi/cli/cmake_variables.py.in" 
               "${CMAKE_CURRENT_SOURCE_DIR}/dune/copasi/cli/cmake_variables.py")

dune_python_install_package(PATH .)