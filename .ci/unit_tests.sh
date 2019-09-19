# unit tests script for all CIs

# make sure we get the right mingw64 version of g++ on appveyor
PATH=/mingw64/bin:$PATH
echo "PATH=$PATH"
echo "MSYSTEM: $MSYSTEM"
echo "DUNECONTROL: ${DUNECONTROL}"
echo "DUNE_OPTIONS_FILE: ${DUNE_OPTIONS_FILE}"
cat ${DUNE_OPTIONS_FILE}

which g++
which python
which cmake
g++ --version
gcc --version
cmake --version

${DUNECONTROL} --opts=${DUNE_OPTIONS_FILE} --only=dune-copasi make --target build_unit_tests
${DUNECONTROL} --opts=${DUNE_OPTIONS_FILE} --only=dune-copasi bexec ctest -j4 -L "unit" --output-on-failure