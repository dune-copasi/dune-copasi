---
id: api
title: Application Programming Interface
sidebar_label: API
---


:::caution Work In Progress
:::
### Doxygen Documentation

The C++ classes are documented in the C++ source code and may be visualized with
[Doxygen](https://www.doxygen.nl/index.html). To locally have the doxygen html
documentation, build `dune-copasi` [manually](install_use.md) and use the
`make doc` target in the build directory. For this, you need to have `doxygen` on your system,

### Preprocessor Definitions

There are two sources of preprocessor definitions:

* **CMake targets**: These are transsitively passed together with the `dune-copas::dune-copasi` target.
  These are stored on the `INTERFACE_COMPILE_DEFINITIONS` property of the target.

  ```bash title="CMakeLists.txt"
  get_target_property(definitions dune-copasi::dune-copasi INTERFACE_COMPILE_DEFINITIONS)
  message("dune-copasi cmake compile definitions: '${definitions}'")
  ```

* **Config header**: We inherit a config header from the dune build system. It's installed under
  `include/dune/copasi/config.h` on the cmake install prefix (e.g. `/opt/dune/`).

The following is a list of useful definitions:

| Definition | Usage |
| ---------| -------------- |
| `DUNE_COPASI_SD_LIBRARY` | Inspect if the precompiled instantiations for single domain models are present. If so, you may skip inclusion of `dune/copasi/model/multidomain_diffusion_reaction.cc` template definitions.
| `DUNE_COPASI_MD_LIBRARY` | Inspect if the precompiled instantiations for multiple domain models are present. If so If so, you may skip inclusion of `dune/copasi/model/diffusion_reaction.cc` template definitions.
| `HAVE_DUNE_COPASI_CONFIG_H` | Whether the `dune/copasi/config.h` header must be included
| `DUNE_USE_FALLBACK_FILESYSTEM` | Whether `dune-copasi` is using [`ghc::filesystem`](https://github.com/gulrak/filesystem)
