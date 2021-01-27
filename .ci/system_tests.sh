#!/usr/bin/env bash

# system tests script for all CIs

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

if test "x$DUNE_OPTIONS_FILE" != "x"; then
  CMAKE_FLAGS="$(. $DUNE_OPTIONS_FILE; eval echo \$CMAKE_FLAGS)"
fi

mkdir dune-copasi/test/build-cmake && cd dune-copasi/test/build-cmake

echo "cmake $CMAKE_FLAGS .."
eval cmake $CMAKE_FLAGS ..

cmake --target build_system_tests
ctest -j4 -L "DUNE_SYSTEMTEST" --output-on-failure
