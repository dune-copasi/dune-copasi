[![Build Status](https://gitlab.dune-project.org/santiago.ospina/dune-copasi/badges/master/pipeline.svg)](https://gitlab.dune-project.org/santiago.ospina/dune-copasi/pipelines)
[![Build Status](https://travis-ci.org/SoilRos/dune-copasi.svg?branch=master)](https://travis-ci.org/SoilRos/dune-copasi)
[![Build status](https://ci.appveyor.com/api/projects/status/6605joy2w17qvca8/branch/master?svg=true)](https://ci.appveyor.com/project/SoilRos/dune-copasi/branch/master)

#### Dependencies

| Software | Version/Branch | Comments |
| ---------| -------------- | -------- |
| muParser | - |
| CMake | 3.10.2 |
| C++ compiler | [C++17](https://en.wikipedia.org/wiki/List_of_compilers#C++_compilers) | 
| [dune-common](https://gitlab.dune-project.org/santiago.ospina/dune-common) | support/dune-copasi
| [dune-logging](https://gitlab.dune-project.org/staging/dune-logging) | master
| [dune-geometry](https://gitlab.dune-project.org/core/dune-geometry) | master
| [dune-grid](https://gitlab.dune-project.org/core/dune-grid) | master
| [dune-uggrid](https://gitlab.dune-project.org/staging/dune-uggrid) | master
| [dune-istl](https://gitlab.dune-project.org/core/dune-istl) | master
| [dune-localfunctions](https://gitlab.dune-project.org/core/dune-localfunctions) | master
| [dune-functions](https://gitlab.dune-project.org/staging/dune-functions) | master
| [dune-typetree](https://gitlab.dune-project.org/santiago.ospina/dune-typetree) | support/dune-copasi
| [dune-pdelab](https://gitlab.dune-project.org/santiago.ospina/dune-pdelab) | support/dune-copasi
| [dune-multidomaingrid](https://gitlab.dune-project.org/santiago.ospina/dune-multidomaingrid) | support/dune-copasi

<!-- 
Preparing the Sources
=========================

Additional to the software mentioned in README you'll need the
following programs installed on your system:

```
  cmake >= 2.8.12
```

Getting started
---------------

If these preliminaries are met, you should run

```
  dunecontrol all
```

which will find all installed dune modules as well as all dune modules
(not installed) which sources reside in a subdirectory of the current
directory. Note that if dune is not installed properly you will either
have to add the directory where the `dunecontrol` script resides (probably
`./dune-common/bin`) to your path or specify the relative path of the script.

Most probably you'll have to provide additional information to `dunecontrol`
(e. g. compilers, configure options) and/or make options.

The most convenient way is to use options files in this case. The files
define four variables:

```
CMAKE_FLAGS      flags passed to cmake (during configure)
```

An example options file might look like this:

```bash
#use this options to configure and make if no other options are given
CMAKE_FLAGS=" \
-DCMAKE_CXX_COMPILER=g++-5 \
-DCMAKE_CXX_FLAGS='-Wall -pedantic' \
-DCMAKE_INSTALL_PREFIX=/install/path" #Force g++-5 and set compiler flags
```

If you save this information into example.opts you can pass the opts file to
dunecontrol via the `--opts option`, e. g.

```bash
  dunecontrol --opts=example.opts all
```

More info
---------

See

```bash
     dunecontrol --help
```

for further options.


The full build system is described in the `dune-common/doc/buildsystem` (Git version) or under `share/doc/dune-common/buildsystem` if you installed DUNE! -->