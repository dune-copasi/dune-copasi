## Testing for dune-copasi

In order to run the tests, you should have built dune-copasi with cmake.

### Unit tests

The unit tests are independent pieces of C++ code tested Google Tests.

```bash
cmake --build . --target build_unit_tests
ctest -L "unit" --output-on-failure
```

### System tests

The system tests are minimal configuration files that are run against the `dune-copasi` executable.

```bash
cmake --build . --target build_system_tests
ctest -L "system" --output-on-failure
```

### Tutorial tests

The tutorial tests are configuration files from the tutorials that are run against the `dune-copasi` executable.

```bash
cmake --build . --target build_tutorial_tests
ctest -L "tutorial" --output-on-failure
```
