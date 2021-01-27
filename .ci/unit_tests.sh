#!/usr/bin/env bash

# unit tests script for all CIs

set -e

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "DUNE_OPTIONS_FILE: ${DUNE_OPTIONS_FILE}"
echo "PWD: $PWD"

which g++
which python
which cmake
g++ --version
gcc --version
cmake --version

if test "x$DUNE_OPTS_FILE" != "x"; then
  CMAKE_FLAGS="$(. $DUNE_OPTS_FILE; eval echo \$CMAKE_FLAGS)"
fi

mkdir dune-copasi/test/build-cmake && cd dune-copasi/test/build-cmake

echo "cmake $CMAKE_FALGS --build .."
eval cmake $CMAKE_FALGS --build ..

cmake --target build_unit_tests
ctest -j4 -L "unit" --output-on-failure
