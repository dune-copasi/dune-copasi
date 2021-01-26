#!/usr/bin/env bash

# build script for all CIs
# requisites
#   * DUNE_OPTIONS_FILE is defined
#   * DUNECONTROL is defined or dunecontrol is in the PATH

set -e

if [[ -z "$DUNECONTROL" ]]
then
  DUNECONTROL="dunecontrol"
fi

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "DUNECONTROL: ${DUNECONTROL}"
echo "DUNE_OPTIONS_FILE: ${DUNE_OPTIONS_FILE}"
cat ${DUNE_OPTIONS_FILE}
echo "PWD: $PWD"

which g++
which python
which cmake
g++ --version
gcc --version
cmake --version

${DUNECONTROL} --opts=${DUNE_OPTIONS_FILE} --only=dune-copasi all
cmake --build dune-copasi/build-cmake --target install
rm -rf dune-copasi/build-cmake
