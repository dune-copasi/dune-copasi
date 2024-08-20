[![Build Status](https://gitlab.dune-project.org/copasi/dune-copasi/badges/master/pipeline.svg)](https://gitlab.dune-project.org/copasi/dune-copasi/pipelines)
[![Build Status](https://github.com/dune-copasi/dune-copasi/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/dune-copasi/dune-copasi/actions/workflows/ci.yml)
[![Netlify Status](https://api.netlify.com/api/v1/badges/6fc6d371-87df-49b5-8e72-e1873fa5d54b/deploy-status)](https://app.netlify.com/sites/dune-copasi/deploys)

<!-- IMPORTANT: Do not add relative links on this file because this is rendered out of source at https://www.dune-project.org/modules/dune-copasi/ too -->

---
# DuneCopasi: A Multi-Compartment Reaction-Diffusion Simulator

DuneCopasi is a C++ library providing a numerical simulator for systems of reaction-diffusion equations on single and multiple compartments using the **D**istributed and **U**nified **N**umerics **E**nvironment [DUNE](https://www.dune-project.org). It supports multi-threaded computations with continuous finite elements on unstructured grids in 1, 2, and 3 dimensions. An application program is made available as a Docker Image, as a Web Assembly binary, or as a GUI through the [Spatial Model Editor](https://spatial-model-editor.github.io) for the manipulation of spatial models of bio-chemical reactions.

# Table of Contents

- [DuneCopasi: A Multi-Compartment Reaction-Diffusion Simulator](#dunecopasi-a-multi-compartment-reaction-diffusion-simulator)
- [Table of Contents](#table-of-contents)
- [Overview](#overview)
  - [Features](#features)
  - [Capabilities](#capabilities)
  - [Model Problem](#model-problem)
- [Getting Started](#getting-started)
  - [Installation](#installation)
  - [Basic Commands](#basic-commands)
  - [Tutorials](#tutorials)
  - [Configuration Options](#configuration-options)
- [Contributing](#contributing)
  - [Reporting Issues](#reporting-issues)
  - [Contributing Guide](#contributing-guide)
  - [Folder Structure](#folder-structure)
- [Additional Information](#additional-information)
  - [License](#license)
  - [Acknowledgements](#acknowledgements)

---

# Overview

DuneCopasi bridges the gap between systems biology and scientific computing by providing a robust tool for simulating complex biochemical networks, especially those with spatial dependencies. It is particularly suitable for researchers and engineers in systems biology, computational biology, and related fields, offering capabilities that include advanced handling of multi-compartment reaction-diffusion systems.

## Features

- **Multi-Compartment Reaction-Diffusion Simulation**: Solve systems of PDEs across multiple compartments with configurable parameters.
- **Continuous Galerkin Method**: Utilize continuous finite elements on unstructured grids in 1D, 2D, and 3D.
- **Multi-Threaded Computation**: Efficient simulation on modern multi-core processors.
- **Interoperability**: Designed to work seamlessly with SBML, TIFF files for initial conditions, and systems biology workflows.
- **Flexible Deployment**: Available as a Docker image, WebAssembly binary, and as a GUI through the Spatial Model Editor.

## Capabilities

- **Configurable Executable**: Single executable with runtime configuration options.
- **SBML Compatibility**: Supports math parsers compatible with the Systems Biology Markup Language (SBML) for model specifications.
- **Custom Grid and TIFF Input**: Accepts custom grids and image files for initial conditions.
- **Advanced Boundary Conditions**: Flexible boundary/transmission conditions for complex biological models.
- **Cross-Diffusion Support**: Allows species to cross-diffuse across compartments.
- **Result Comparison**: Built-in tools for comparing simulation results with objective functions.
- **Random Field Generator**: Generates statistical spatial variations using embedded random field algorithms.

## Model Problem

DuneCopasi solves the following system of partial differential equations (PDEs):

$$
\begin{aligned}
  \partial_t (\phi_{ik} u_{ik}) &= \nabla \cdot \vec{\mathcal{j}}_{ik}+\mathcal{R}_{ik}(\mathbf{u})\qquad &\text{in }\Omega_k\subset \mathbb{R}^d,\\
  u_{ik} &= u^{(0)}_{ik} \qquad &\text{on }\Gamma_k^D\subseteq\partial\Omega_k,\\
  \vec{\mathcal{j}}_{ik}\cdot \mathbf{n}_k &= \sum_{l\in T_k}\mathcal{T}_{ikl}(\mathbf{u}) \qquad &\text{on }\partial\Omega_k \setminus \Gamma^D_{k},\\
  \vec{\mathcal{j}}_{ik} &:= \sum_{j\in N_k}\mathsf{D}_{ijk}\nabla u_{jk}\\
\end{aligned}
$$

where $i\in N_k$ is the subscript for each mass balance equation in the $k$-th compartment and with $k\in K$, $T_k:=\{l\in K : \Gamma_{kl}\neq\varnothing\}$, $\mathbf{u}_k:=(u_{1k},\ldots,u_{\text{dim}(N_k)k})$, $\mathbf{u} := (\mathbf{u}_1, \ldots, \mathbf{u}_{\text{dim}(K)})$, $\mathbf{n}_k$ the unit outer normal on $\Omega_k$, $\overline{\Omega}:=\bigcup_{k\in K}\overline{\Omega}_k$, and $\Gamma_{kl}$ is $\partial\Omega_k\cap\overline{\Omega}_l$ if $k \neq l$ or $\partial\Omega_k\cap\overline{\Omega}$ otherwise.


All of the parameterization for the equations are easily configurable at run-time using the command line or a configuration file. That is, the underlying grid, the partition for each sub-domain $\Omega_k$ on the grid, the reaction operator $\mathcal{R}_{ik}(\mathbf{u})$, the non-linear transmission conditions $\mathcal{T}_{ikl}(\mathbf{u})$ on $\Gamma_{kl}$, the Dirichlet boundary conditions $u^{(0)}_{ik}$ on $\Gamma^D_{k}$, the storage terms $\phi_{ik}$, and the cross-diffusion terms $\mathsf{D}_{ijk}$.

# Getting Started

## Installation

For detailed installation instructions specific to each version, please refer to our comprehensive [online documentation](https://dune-copasi.netlify.app).

## Basic Commands

```
USAGE: dune-copasi [options]

Numerical simulator for diffusion-reaction systems on single or multiple compartments

Website: <https://dune-copasi.netlify.app>

OPTIONS:

    -h / --help          - Display this help.
    --help-full          - Display this help with long descriptions.
    --version            - Display the version of this program.
    --parser-default     - Display the default parser.
    --parser-list        - Display the parsers available for this program.
    --dimension-list     - Display the grid dimensions available for this program.
    --dump-config        - Dumps configuration in the INI format to stdout.
    --config=<string>    - Specifies a config file in INI format. See Configuration Options.
    --{key}={value}      - Overrides key=value sections of the config file. See Configuration Options.
```

## Tutorials

To help you get started with DuneCopasi, we offer a dedicated [tutorial website](https://dune-copasi.netlify.app/tutorials) featuring step-by-step guides and practical examples. These tutorials are designed to walk you through various aspects of the software.

## Configuration Options

We have a rich set of options to configure many possible settings for our problem model.
You can find detailed information of each option by invoking the `--help-full` command or by exploring our [online documentation](https://dune-copasi.netlify.app).

# Contributing

## Reporting Issues

If you encounter any issues, please report them on the [issue tracker](https://gitlab.dune-project.org/copasi/dune-copasi/issues) or by contacting us via the email ![Service Desk](https://img.shields.io/badge/dune--copasi%40dune--project.org-informational).

## Contributing Guide

We welcome contributions from the community. Please read our Contributing Guide in `CONTRIBUTING.md` file for details on our code of conduct, and the process for submitting merge requests.

## Folder Structure

The following may help you understand the general structure of this project.
Below, you can find the folder structure with the most prominent parts of this repository and a short description of their purpose.
Notice that many of these sub-directories have their own `README.md` file explaining how to use its contents.

```
.                                          # . Root directory
├── .ci                                    # ├── Scripts to run during Continuous Integration (CI)
├── .github                                # ├── GitHub configuration files
├── .gitlab                                # ├── GitLab configuration files
├── cmake                                  # ├── CMake scripts and modules
├── doc                                    # ├── Documentation files
│   ├── docusaurus                         # │   ├── Files for online documentation
│   │   ├── docs                           # │   │   ├── Main online documentation files
│   │   │   ├── assets/ini                 # │   │   │   ├── Configuration examples for online documentation (tested)
│   │   │   └── assets/config_opts.json    # │   │   │   └── Configuration options database (source of truth for all docs)
│   │   ├── tutorials                      # │   │   ├── Main online tutorial files
│   │   │   └── assets/ini                 # │   │   │   └── Configuration examples for online tutorials (tested)
│   │   └── ...                            # │   │   └── [...]
│   ├── doxygen                            # │   ├── Configuration files for autogenerated C++ documentation
│   ├── joss                               # │   ├── Paper submitted to the Journal of Open Source Software
│   └── man                                # │   └── Documentation for `man`
├── docker                                 # ├── Docker files for different environments
├── dune/copasi                            # ├── C++ header files
│   ├── common                             # │   ├── Common infrastructure
│   ├── concepts                           # │   ├── Type constrains/requirements
│   ├── finite_element                     # │   ├── Finite element utilities
│   ├── grid                               # │   ├── Grids utilities
│   ├── model                              # │   ├── Simulation logic for a model problem
│   │   └── diffusion_reaction             # │   │   └── Diffusion-Reaction specifics of the model problem
│   ├── parser                             # │   ├── Interface and implementation of math parsers
│   └── solver/istl                        # │   └── Custom linear algebra solvers/preconditioners for dune-istl
├── npm                                    # ├── WebAssembly package files for NPM
├── src                                    # ├── Object files and main program (mirrors dune/copasi folder)
├── test                                   # ├── Testing suite
│   ├── data                               # │   ├── Custom data used for tests
│   ├── dune                               # │   ├── Unit tests on specific C++ constructs (mirrors dune/copasi folder)
│   └── *.ini                              # │   └── Configuration files for system tests
├── util                                   # ├── Misc utilities to handle the repository
├── .gitlab-ci.yml                         # ├── GitLab CI configuration file
├── CHANGELOG.md                           # ├── List of major changes
├── CMakeFiles.txt                         # ├── Main CMake configuration file
├── CODE_OF_CONDUCT.md                     # ├── Code of conduct
├── CONTRIBUTING.md                        # ├── Contribution guilde
├── Dockerfile                             # ├── Main docker file to build the project
├── dune-copasi.opts                       # ├── DUNE options file with default CMake options
├── LICENSE.md                             # ├── BSD 2-Clause License file
├── README.md                              # ├── Overview documentation of repository
└── ...                                    # └── [...]
```

# Additional Information

## License

This project is licensed under the BSD 2-Clause License - see the `LICENSE.md` file for details.

## Acknowledgements

Research funded by the German federal Ministry of Education and Research (BMBF) FKZ 031L0158
