#!/usr/bin/env bash

# system tests script for all CIs

set -e

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

mkdir dune-copasi/test/build-cmake && cd dune-copasi/test/build-cmake
cmake --build ..
cmake --target build_system_tests
ctest -j4 -L "DUNE_SYSTEMTEST" --output-on-failure
