---
id: install_use
title: Installation and Usage
sidebar_label: Install & Use
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

:::caution Work In Progres
:::

There are two forms to use `DuneCopasi`. The first and easiest one is using a
[Docker Container](https://www.docker.com/) where the software is boundled such
that no installation other than docker is required. On the second one, the
compilation of the source code is compulsory.

## Containerized Usage

In this mode, there is no need to fullfil specific requirements other than
those for docker installation.

### Install Docker
First, get and install Docker following the
[docker installation instructions](https://docs.docker.com/get-docker/).

### Pull Container
Then, you may get the latest stable `DuneCopasi` container by pulling it from
the [GitLab registry](https://gitlab.dune-project.org/copasi/dune-copasi/container_registry)

```bash
docker pull registry.dune-project.org/copasi/dune-copasi/dune-copasi:latest
```

### Run `DuneCopasi`
Finally, run the `dune_copasi_md` executable from the container
```bash
# TODO
```

## Manual Installation

Locally building and installing `DuneCopasi` requires to first obtain the
dependencies, then compile the dependent dune modules and finally installing on
a final directory.

### Dependencies

The following list of software is required to install and use `DuneCopasi`:

:::info
Notice that some required dune modules are forks of original reopsitories and
are placed under the [COPASI namespace](https://gitlab.dune-project.org/copasi/)
on the [DUNE GitLab](https://gitlab.dune-project.org/).
:::

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
| [COPASI/dune-logging](https://gitlab.dune-project.org/copasi/dune-logging)                  | support/dune-copasi |
| [COPASI/dune-typetree](https://gitlab.dune-project.org/copasi/dune-typetree)                | support/dune-copasi |
| [COPASI/dune-pdelab](https://gitlab.dune-project.org/copasi/dune-pdelab)                    | support/dune-copasi |
| [COPASI/dune-multidomaingrid](https://gitlab.dune-project.org/copasi/dune-multidomaingrid)  | support/dune-copasi |

### Installation

The first four can be obtained by your prefered package manager in unix-like operating systems. e.g.

<Tabs
  groupId="operating-systems"
  defaultValue="win"
  values={[
      {label: 'Debian/Ubuntu', value: 'deb', },
      {label: 'macOS', value: 'mac', },
    ]
  }>

  <TabItem value="deb">

```bash
apt update
apt install cmake gcc g++ libtiff-dev libmuparser-dev git
```

  </TabItem>
  <TabItem value="mac">

```bash
# Apple Command Line Tools and Brew assumed to be available
brew update
brew install cmake gcc libtiff muparser git
```

  </TabItem>
</Tabs>

The required `DUNE` modules (including `DuneCopasi`) can be obtained via internet
by using [`git`](https://git-scm.com/). For smooth installation, is better place
all the dune modules within the same directory.

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
git clone -b support/dune-copasi --recursive https://gitlab.dune-project.org/copasi/dune-logging
git clone -b support/dune-copasi https://gitlab.dune-project.org/copasi/dune-typetree
git clone -b support/dune-copasi https://gitlab.dune-project.org/copasi/dune-pdelab
git clone -b support/dune-copasi https://gitlab.dune-project.org/copasi/dune-multidomaingrid
git clone -b v0.2.0 https://gitlab.dune-project.org/copasi/dune-copasi
```

Then build and install the `DUNE` modules with the `dunecontrol` script:
```bash
# configure and build dune modules
./dune-common/bin/dunecontrol make all

# install dune-copasi (this operation may requiere sudo)
./dune-common/bin/dunecontrol bexec make install

# remove source and build files
cd .. && rm -r ~/dune-modules
```

For further info on dune module installation process, please check out
the [dune-project web page](https://www.dune-project.org/doc/installation/).

:::tip Custom installation path
If you don't want to install the `DUNE` modules system-wide, you can set the
`CMAKE_INSTALL_PREFIX` to a custom directory. See CMake documentation
[here](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html).
:::

### Run `DuneCopasi`

In case you installed `DuneCopasi` system-wide, you should be able to call the
program `dune_copasi_md` from your command line accompained with a configuration
file. Otherwise, the executable will be under `bin/` folder on the installed
directory.

```bash
# TODO
```

To find out the appropiated contents on the configuration file, check out
the [Parameter Tree](param_tree.md) documentation.
