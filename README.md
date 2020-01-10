[![Build Status](https://gitlab.dune-project.org/copasi/dune-copasi/badges/master/pipeline.svg)](https://gitlab.dune-project.org/copasi/dune-copasi/pipelines)
[![Build Status](https://travis-ci.org/SoilRos/dune-copasi.svg?branch=master)](https://travis-ci.org/SoilRos/dune-copasi)
[![Build status](https://ci.appveyor.com/api/projects/status/6605joy2w17qvca8/branch/master?svg=true)](https://ci.appveyor.com/project/SoilRos/dune-copasi/history)

# dune-copasi

Solver for reaction-diffusion systems in multiple compartments

 * Solve a reaction-diffusion system for each comartment
 * Each compartment may have different system with different number of variables
 * Neumann flux at the interface of compartments for variables with
   the same name on the two compartments
 * Easy to modify configuration file
 * Initial conditions might be either TIFF file or a math expression.
 * Solved using the finite element method
 * Output in the VTK format
 * Currently it only supports 2D simulations

This project is made under the umbrella of the 
[*Distributed and Unified Numerics Environment* `DUNE`](https://www.dune-project.org/) and the
[*Biochemical System Simulator* `COPASI`](http://copasi.org/). 
Altought the rationale of the design is always driven by biochemical process (e.g. cell biology), 
this software is not limited to this scope and can be used for other processes involving reaction-diffusion systems.

## Graphical User Interface for SMBL files

For those working in bio-informatics there exist a grafical user interface for SMBL files!
The GUI is able to convert non-spatial SBML models of bio-chemical reactions into 
2d spatial models, and to simulate them with `dune-copasi`.

https://github.com/lkeegan/spatial-model-editor

## Installation

This requires that you have installed the following packages before the actual installation of `dune-copasi`.

| Software | Version/Branch | Comments |
| ---------| -------------- | -------- |
| [CMake](https://cmake.org/) | 3.1 |
| C++ compiler | [C++17](https://en.wikipedia.org/wiki/List_of_compilers#C++_compilers) | 
| [libTIFF](http://www.libtiff.org/) | 3.6.1 |
| [muParser](https://beltoforion.de/article.php?a=muparser) | 2.2.5 |
| [dune-common](https://gitlab.dune-project.org/copasi/dune-common) | releases/2.7 | https://gitlab.dune-project.org/core/dune-common
| [dune-geometry](https://gitlab.dune-project.org/core/dune-geometry) | releases/2.7 | https://gitlab.dune-project.org/core/dune-geometry
| [dune-grid](https://gitlab.dune-project.org/core/dune-grid) | releases/2.7 | https://gitlab.dune-project.org/core/dune-grid
| [dune-uggrid](https://gitlab.dune-project.org/staging/dune-uggrid) | releases/2.7 | https://gitlab.dune-project.org/staging/dune-uggrid
| [dune-istl](https://gitlab.dune-project.org/core/dune-istl) | releases/2.7 | https://gitlab.dune-project.org/core/dune-istl
| [dune-localfunctions](https://gitlab.dune-project.org/core/dune-localfunctions) | releases/2.7 | https://gitlab.dune-project.org/core/dune-localfunctions
| [dune-functions](https://gitlab.dune-project.org/staging/dune-functions) | releases/2.7 | https://gitlab.dune-project.org/staging/dune-functions
| [dune-logging](https://gitlab.dune-project.org/staging/dune-logging) | support/dune-copasi (recursive) | https://gitlab.dune-project.org/copasi/dune-logging
| [dune-typetree](https://gitlab.dune-project.org/copasi/dune-typetree) | support/dune-copasi | https://gitlab.dune-project.org/copasi/dune-typetree
| [dune-pdelab](https://gitlab.dune-project.org/copasi/dune-pdelab) | support/dune-copasi | https://gitlab.dune-project.org/copasi/dune-pdelab
| [dune-multidomaingrid](https://gitlab.dune-project.org/copasi/dune-multidomaingrid) | support/dune-copasi | https://gitlab.dune-project.org/copasi/dune-multidomaingrid

The first four can be obtained by your prefered package manager in unix-like operating systems. e.g.

```bash
# if you are in a debian/ubuntu OS
apt update
apt install cmake gcc g++ libtiff-dev libmuparser-dev git

# if you are in a macOS
xcode-select --install # Apple Command Line Tools
brew update
brew install cmake gcc libtiff muparser git
```

Now, the dune modules (including `dune-copasi`) can be all checkout in a same folder and be installed in one go. 

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
git clone -b master https://gitlab.dune-project.org/copasi/dune-copasi

# configure and build dune modules
./dune-common/bin/dunecontrol make all

# install dune-copasi (this operation may requiere sudo)
./dune-common/bin/dunecontrol --only=dune-copasi bexec make install

# if you do not want to install dune-copasi system-wide, you can set
# the CMAKE_INSTALL_PREFIX to a non restricted folder
# see https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html

# remove source and build files
cd .. && rm -r ~/dune-modules
```

For further info on dune module installation process, please check out 
the [dune-project web page](https://www.dune-project.org/doc/installation/)

## Usage 

If you installed `dune-copasi` system-wide, you should be able to call the program
`dune_copasi` from your command line accompained with a configuration file.

```bash
dune_copasi config.ini
```

### Configuration File

The configuration file follows [INI file format](https://en.wikipedia.org/wiki/INI_file).
It should contain at least two sections: `grid`, and `model`, whereas a third section 
`logging` is optional.

#### Grid

The grid section is fairly simple as it only contains the path to a [gmsh file](http://gmsh.info/)
and its initial refinement level:

```ini
[grid]
file = my_gmsh_file.msh
initial_level = 1
```

#### Model

The model section starts with the definitions of the simulation time interval 
and the polynomal order of the finite element method (currently only supports `1` and `2`).

```ini
[model]
begin_time = 0.
end_time = 10
time_step = 0.1
order = 1
```

The following is the definition of the compartments of the model. 
Each compartment corresponds to a *physical group* in the gmsh jargon.
Although the gmsh format allows you to name such physical groups, 
we still need to assign them to a `dune-copasi` compartmet and for that
we use the *physical group* index. Notice that uses 0-based indices while 
`gmsh` uses 1-based indices. In other words, 
`<gmsh_physical_group> = <dune_copasi_compartment> - 1`. 
Let's say for example that there is two *physical groups* in our gmsh file
and we are going to use them as `nucleous` and `cytoplasm` compartments:

```ini
[model.compartments]
 # nucleous corresponds to the physical group 1 in the gmsh file
nucleous  = 0
 # cytoplasm corresponds to the physical group 2 in the gmsh file
cytoplasm = 0
```

Now, each of these compartments will define its own initial conditions,
its diffusion-reaction system, and its vtk writter. For that, you have to expand
the `model` section with the defined compartments, e.g. `model.nucleous` or `model.cytoplasm`.
The subsection `initial`, `reaction`, `diffusion` and `operator` define the system variables 
and its properties. You can put as many variables as desired as long as they are the same 
in this three subsections. Notice that each variable defines a new diffusion-reaction partial 
differential equation associated with it. 

  * The `initial` subsection allow the initialization of each of the variables. 
  * The `diffusion` subsection defines the diffusion coefficient math 
  expression associated with each variable. It may only depend on the grid coordinates `x` and `y`.
  * The `reaction` subsection defines the reaction part of each equation in the PDE. 
  Since this is the souce of non-linearities, it allows to be dependent on other defined variables
  within the compartment. This section has to include yet another subsection called `jacobian`.
  * The `reaction.jacobian` subsection must describe the 
  [jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
  of the `reaction` part. It must follow the syntax of `d<var_i>_d<var_j>`, which 
  reads as the *partial derivate of `<var_i>` with respect to `<var_j>`*. 
  * The `operator` subsection is an experimental feature and we recommend to set all variables to the 
  same index, e.g. 0. 
  * Finally, the subsection `writer` will define the file name for the vtk output. 
  
For example, the following `mode.nucleous` section defines a [Gray-Scott
model with `F=0.0420` and `k=0.0610`](http://mrob.com/pub/comp/xmorphia/F420/F420-k610.html):

```ini
[model.nucleous.initial]
u = 0.5
v = (x>0) && (x<0.5) && (y>0.) && (y<0.5) ? 1 : 0

[model.nucleous.diffusion]
u = 2e-5
v = 2e-5/2

[model.nucleous.reaction]
u = -u*v^2+0.0420*(1-u)
v = u*v^2-(0.0420+0.0610)*v

[model.nucleous.reaction.jacobian]
du_du = -v^2-0.0420
du_dv = -2*u*v
dv_du = v^2
dv_dv = 2*u*v-(0.0420+0.0610)

[model.nucleous.operator]
u = 0
v = 0

[model.nucleous.writer]
file_name = nucleous_output
```
The `model.cytoplasm` would have to be defined in similar way. An important aspect when working 
with different compartments is the interface fluxes. In `dune-copasi` thex fluxes are set 
automatically to [dirichlet-dirichlet](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition)
boundary conditions iff the variable is shared between the two intersecting domains. Further 
improvements will come for interface treatment.

#### Logging

The logging settings are directly forwarded to the `dune-logging` module. Please check 
its doxygen documentation for detailed information. A simple configuration is the following:

```ini
[logging]
# possible levels: off, critical, error, waring, notice, info, debug, trace, all
default.level = trace

[logging.sinks.stdout]
pattern = [{reldays:0>2}-{reltime:8%T}] [{backend}] {msg}
```