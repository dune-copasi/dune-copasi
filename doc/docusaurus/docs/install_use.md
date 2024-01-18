---
id: install_use
title: Installation and Usage
sidebar_label: Install & Use
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

## Usage

There are three different ways to use `dune-copasi`:

* [GUI: Graphical User Interface](#graphical-user-interface)
* [INI: Configuration File](#configuration-file)
* [API: Application Programming Interface](#application-programming-interface)

### Graphical User Interface

The [Spatial Model Editor](https://github.com/spatial-model-editor/spatial-model-editor)
is a **user friendly** GUI editor to create and edit 2D
spatial *Systems Biology Markup Language* ([SBML](https://en.wikipedia.org/wiki/SBML))
models of bio-chemical reactions. Additionally, it can simulate them with `dune-copasi`. A big
advantage of this package is that is tailored for biologists and is available
with just a pair of clicks on the major plataforms. Find more information
[here](https://spatial-model-editor.readthedocs.io/en/latest/quickstart/get-started.html)!

[![sme](./img/spatial-model-editor.png)](https://spatial-model-editor.readthedocs.io/en/latest/quickstart/get-started.html)

### Configuration File

In this form, `dune-copasi` provides one executable for single compartment
systems (`dune-copasi-sd`) and another one for multiple compartment systems
(`dune-copasi-md`). Both executables expect one [`INI` configuration file](ini_file.md)
which shall contain all the information to perform the simulation.

```bash
dune-copasi-md config.ini
```

Find more information about available configurations on the
[Parameter Tree](param_tree.md) documentation. This form may be
installed in one of the following procedures:

* [Docker Runner](#docker-runner)
* [Debian/Ubuntu and macOS Packages](#debianubuntu-and-macos-packages)
* [Manual Installation](#manual-installation)

### Application Programming Interface

The `dune-copasi` C++ objects may be consumed by other programs in order
to generate custom simulation rules, to couple intermediate steps with other
tools, or to implement another GUI, etc. In such a case, `dune-copasi` must be
available in development mode and the downstream library is expected to
[consume the library](#importing-cmake-targets) by using the
[CMake build system](https://cmake.org) and use the [C++ objects](api.md) in code.
This form is available on:

* [Manual Installation](#manual-installation)
* [Docker Build](#docker-build)

## Installation

### Graphical User Interface

To install the Spatial Model Editor, please refer to its installation
instructions:

* https://spatial-model-editor.readthedocs.io/en/stable/quickstart/get-started.html#installation

### Docker Runner

The easiest form to use our executables for INI usage, is by using a
[Docker Container](https://www.docker.com/). There, the software is boundled
such that no installation other than docker is required.

#### Install Docker
First, get and install Docker following the
[docker installation instructions](https://docs.docker.com/get-docker/).

#### Prepare a working directory

To be able to easily share data between your operating system and the docker
container, prepare a working directory with read/write rights to _other users_
(e.g. a folder named `dune-copasi`):

```bash
mkdir -m o+rw dune-copasi && cd dune-copasi
```

This working directory will be accessible to your text editor, paraview as
well as to the `dune-copasi-md` executable inside the docker container. Thus,
move or create your configuration files into it at will.

```
# setup/write ini and input files for dune-copasi...
nano config.ini
```

#### Run the program

Here, you may pull and run the latest stable container from our
[GitLab registry](https://gitlab.dune-project.org/copasi/dune-copasi/container_registry).
To do so, call the docker container with a configuration
file `config.ini` using one of the following commands on the terminal:

<Tabs
  groupId="docker-run"
  defaultValue="md"
  values={[
      {label: 'Multi Domain', value: 'md', },
      {label: 'Single Domain', value: 'sd', },
    ]
  }>

  <TabItem value="md">

```bash
docker run -v $PWD \
  registry.dune-project.org/copasi/dune-copasi/dune-copasi:latest\
  config.ini
```

  </TabItem>
  <TabItem value="sd">

```bash
docker run -v $PWD \
  --entrypoint=dune-copasi-sd \
  registry.dune-project.org/copasi/dune-copasi/dune-copasi:latest\
  config.ini
```

  </TabItem>
</Tabs>

The results of those computations will be written on current
directory as mentioned above. For more information about running docker images,
visit the [`docker run` documentation](https://docs.docker.com/engine/reference/run/).

### Debian/Ubuntu and macOS Packages

For Debian/Ubuntu/macOS users that want to make use of `dune-copasi`
with [INI](#configuration-file) usage, installation is as simple as:

<Tabs
  groupId="package-manager"
  defaultValue="apt"
  values={[
      {label: 'Debian/Ubuntu (apt)', value: 'apt', },
      {label: 'macOS/Linux (brew)', value: 'brew', },
    ]
  }>

  <TabItem value="apt">

```bash
curl -fsSL https://gitlab.dune-project.org/copasi/dune-copasi/-/jobs/artifacts/latest/raw/packages/dune-copasi-runtime.deb?job=build:debian_clang -o dune-copasi-runtime.deb
apt install ./dune-copasi-runtime.deb
```

  </TabItem>
  <TabItem value="brew">

```bash
brew tap dune-copasi/tap
brew install dune-copasi
```

  </TabItem>
</Tabs>


Once installed, the programs `dune-copasi-sd` and `dune-copasi-md` will be
available on the command line:

```bash
dune-copasi-md config.ini
```

To remove the package call the following command on the terminal

<Tabs
  groupId="package-manager"
  defaultValue="apt"
  values={[
      {label: 'Debian/Ubuntu (apt)', value: 'apt', },
      {label: 'macOS/Linux (brew)', value: 'brew', },
    ]
  }>

  <TabItem value="apt">

```bash
apt remove dune-copasi-runtime
```

  </TabItem>
  <TabItem value="brew">

```bash
brew uninstall dune-copasi
```

  </TabItem>
</Tabs>

### Docker Build

Advanced users, who may want to make modifications the `dune-copasi`
code but do not to install all the dependencies may opt for a docker build. In
this case, you must download the `dune-copasi` source code, modify it, and build
a new local docker image:

```bash
# fetch source code from git
git clone https://gitlab.dune-project.org/copasi/dune-copasi

# enter dune-copasi directory
cd dune-copasi

# checkout the branch you want to modify (e.g. latest)
git checkout latest

# modify source code at will
# ...

# build a new docker image from modified code (tag: dune-copasi)
docker build -t dune-copasi .
```

This will build all dune dependencies as well as the new modified version of
`dune-copasi`. Then, follow the [Docker runner guide](#docker-runner) to run
the new image with tag `dune-copasi`.

### Manual Installation

Finally, to locally build and install `dune-copasi` we requires to obtain, compile
and install a variety of dependencies.

#### Operating Systems and Compilers

The manual installation is known to compile and run under Debian/Ubuntu, macOS, and Windows.
In all three cases with [Clang](https://clang.llvm.org) and [GCC](https://gcc.gnu.org)
based compilers. It is however, known to not be compatible with Microsoft Visual Studio.


#### Dependencies

The following list of software is required to install and use `dune-copasi`:

| Software | Version/Branch |
| ---------| -------------- |
| [CMake](https://cmake.org/)                                                                 | >= 3.1 |
| C++ compiler  | >= [C++17](https://en.wikipedia.org/wiki/List_of_compilers#C++_compilers) |
| [libTIFF](http://www.libtiff.org/)                                                          | >= 3.6.1 |
| [muParser](https://beltoforion.de/article.php?a=muparser)                                   | >= 2.2.5 |
| [dune-common](https://gitlab.dune-project.org/copasi/dune-common)                           | == 2.7 |
| [dune-geometry](https://gitlab.dune-project.org/core/dune-geometry)                         | == 2.7 |
| [dune-grid](https://gitlab.dune-project.org/core/dune-grid)                                 | == 2.7 |
| [dune-uggrid](https://gitlab.dune-project.org/staging/dune-uggrid)                          | == 2.7 |
| [dune-istl](https://gitlab.dune-project.org/core/dune-istl)                                 | == 2.7 |
| [dune-localfunctions](https://gitlab.dune-project.org/core/dune-localfunctions)             | == 2.7 |
| [dune-functions](https://gitlab.dune-project.org/staging/dune-functions)                    | == 2.7 |
| [dune-logging](https://gitlab.dune-project.org/staging/dune-logging)                        | == 2.7 |
| [dune-multidomaingrid](https://gitlab.dune-project.org/extensions/dune-multidomaingrid)     | == 2.7 |
| [COPASI/dune-typetree](https://gitlab.dune-project.org/copasi/dune-typetree)                | `support/dune-copasi-latest` |
| [COPASI/dune-pdelab](https://gitlab.dune-project.org/copasi/dune-pdelab)                    | `support/dune-copasi-latest` |

:::info
Notice that some required dune modules are forks of original repositories and
are placed under the [COPASI namespace](https://gitlab.dune-project.org/copasi/)
on the [DUNE GitLab](https://gitlab.dune-project.org/).
:::

#### Dune Options File

An important part of the installation procedure is to tune the build system
flags to accommodate the build to your system. This is done via the [dune options
file](https://dune-project.org/doc/installation/#storing-flags-in-an-options-file).
In essence, is just a bash script that sets different flags (mainly the flags
for CMake `CMAKE_FLAGS`). While the dune project usually leaves this open for the user, we
provide a `dune-copasi.opts` file with sensible default options for the main
operating systems.

This file can be called to show the current configuration:

```bash
./dune-copasi.opts
```

There are two form to add flags to the `dune-copasi.opts` file:

1. Setting environmental variables starging with `CMAKE_` and `DUNE_`.
   These variables will be automatically included into the list of cmake flags.

  ```bash
  export CMAKE_INSTALL_PREFIX="/path/to/install/prefix/"
  export CMAKE_GENERATOR="'Ninja'"
  ./dune-copasi.opts
  ```

2. Appending new flags into the `CMAKE_FLAGS` variable in the `dune-copasi.opts` file
  (care must be taken to scape quotes):

  ```bash
  echo 'CMAKE_FLAGS+=" -DCMAKE_INSTALL_PREFIX=/path/to/install/prefix/"' >> dune-copasi.opts
  echo "CMAKE_FLAGS+=\" -DCMAKE_GENERATOR='Ninja'\"" >> dune-copasi.opts
  ./dune-copasi.opts
  ```

For more information about the possible options and the dune options file, check out
the dune [installation documentation](https://dune-project.org/doc/installation/)
and the [build system documentation](https://dune-project.org/sphinx/core-2.7/).

#### Installation

The first four dependencies can be obtained by your preferred package manager in
unix-like operating systems. e.g.

<Tabs
  groupId="package-manager"
  defaultValue="apt"
  values={[
      {label: 'Debian/Ubuntu (apt)', value: 'apt', },
      {label: 'macOS/Linux (brew)', value: 'brew', },
    ]
  }>

  <TabItem value="apt">

```bash
apt update
apt install cmake gcc g++ libtiff-dev libmuparser-dev git
```

  </TabItem>
  <TabItem value="brew">

```bash
brew update
brew install cmake gcc libtiff muparser git
```

  </TabItem>
</Tabs>

The required `DUNE` modules (including `dune-copasi`) can be obtained via
internet by using [`git`](https://git-scm.com/). For smooth installation, is
better place all the dune modules within the same directory.

:::caution DUNE Version incompatibility
Notice that this procedure assumest that you don't have previous installation of DUNE in
your system. If that's the case, uninstall DUNE before continouing, otherwise, version
conflicts may be difficult if not impossible to resolve.
:::

```bash
# prepare a folder to download and build dune modules
mkdir ~/dune-modules && cd ~/dune-modules

# fetch dependencies & dune-copasi in ~/dune-modules folder
git clone -b releases/2.7 https://gitlab.dune-project.org/core/dune-common
git clone -b releases/2.7 https://gitlab.dune-project.org/core/dune-geometry
git clone -b releases/2.7 https://gitlab.dune-project.org/core/dune-grid
git clone -b releases/2.7 https://gitlab.dune-project.org/staging/dune-uggrid
git clone -b releases/2.7 https://gitlab.dune-project.org/core/dune-istl
git clone -b releases/2.7 https://gitlab.dune-project.org/core/dune-localfunctions
git clone -b releases/2.7 https://gitlab.dune-project.org/staging/dune-functions
git clone -b releases/2.7 https://gitlab.dune-project.org/extensions/dune-multidomaingrid
git clone -b releases/2.7 --recursive https://gitlab.dune-project.org/staging/dune-logging
git clone -b support/dune-copasi-latest https://gitlab.dune-project.org/copasi/dune-typetree
git clone -b support/dune-copasi-latest https://gitlab.dune-project.org/copasi/dune-pdelab
git clone -b latest https://gitlab.dune-project.org/copasi/dune-copasi

# apply patches
git apply -C dune-common dune-copasi/.ci/dune-common.patch
```

Then, build and install the `DUNE` modules with the `dunecontrol` script:

```bash
# choose an installation path      (this is read by the 'dune-copasi.opts' file)
export CMAKE_INSTALL_PREFIX=/opt/dune/

# configure and build dune modules                                (go grab a coffee)
./dune-common/bin/dunecontrol --opts=dune-copasi/dune-copasi.opts all

# install binaries and libraries                                  (may require sudo)
./dune-common/bin/dunecontrol --opts=dune-copasi/dune-copasi.opts bexec make install

# remove source and build files
cd ~ && rm -r ~/dune-modules

# include dune binaries into your path
echo "export PATH=${CMAKE_INSTALL_PREFIX}/bin:\$PATH" >> $HOME/.bashrc
```

For further information on dune module installation process, please check out
the [dune-project web page](https://www.dune-project.org/doc/installation/).

#### Run the program

Once installed, the programs `dune-copasi-sd` and `dune-copasi-md` will be
available on the command line:

```bash
dune-copasi-md config.ini
```
#### Importing CMake targets

If you additionally want to use the [API](#application-programming-interface) for
development, you must find and consume the CMake targets from `dune-copasi`
in your project as follows:

```
# ...
find_package(dune-copasi IMPORTED REQUIRED)
target_link_libraries(my_app PRIVATE dune-copasi::dune-copasi)
# ...
```

If `dune-copasi` was installed on a custom directory
(e.g. using `CMAKE_INSTALL_PREFIX=/opt/dune`), it may be possible that you need to
pass such directory to the `CMAKE_PREFIX_PATH` when building the project. This way,
CMake can find our targets and configuration:

```bash
cmake -DCMAKE_PREFIX_PATH:PATH=/opt/dune /path/to/app/source/
```
