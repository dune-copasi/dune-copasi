# unit tests script for all CIs

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "PWD: $PWD"

which g++
which python
which cmake
g++ --version
gcc --version
cmake --version

mkdir dune-copasi/test/cmake-build && cd dune-copasi/test/cmake-build
cmake --build ..
cmake --target build_unit_tests
ctest -j4 -L "unit" --output-on-failure
