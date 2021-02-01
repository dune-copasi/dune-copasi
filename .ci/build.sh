#!/usr/bin/env bash

# build script for all CIs
# requisites
#   * DUNE_OPTIONS_FILE is defined

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "DUNE_OPTIONS_FILE: ${DUNE_OPTIONS_FILE}"
cat ${DUNE_OPTIONS_FILE}
echo "PWD: $PWD"

which g++
which python
which cmake
g++ --version
gcc --version
cmake --version

CMAKE_FLAGS="$(. ${DUNE_OPTIONS_FILE}; eval echo \$CMAKE_FLAGS)"

mkdir dune-copasi/build-cmake && cd dune-copasi/build-cmake

echo "cmake $CMAKE_FLAGS .."
eval cmake $CMAKE_FLAGS ..

cmake --build . --target install
