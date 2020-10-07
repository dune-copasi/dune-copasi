---
id: changelog
title: Changelog
sidebar_label: Changelog
---

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--
Guiding Principles

    Changelogs are for humans, not machines.
    There should be an entry for every single version.
    The same types of changes should be grouped.
    Versions and sections should be linkable.
    The latest version comes first.
    The release date of each version is displayed.
    Mention whether you follow Semantic Versioning.

Types of changes

    Added         for new features.
    Changed       for changes in existing functionality.
    Deprecated    for soon-to-be removed features.
    Removed       for now removed features.
    Fixed         for any bug fixes.
    Security      in case of vulnerabilities.
 -->

## [Unreleased]
### Added
- User documentation !28
- Scape VTK write by omitting `[writer]` section !24
- Added Explicit and Diagonally Runge-Kutta time solvers !25
- 3D simulations are now possible !27
### Changed
- Improved errors when a config file is incorrect !25
- Update meaningful logging values on the test ini file !25
- Clean up logger output !25
- Time stepping is now done outside the model class !25
- Writer might write on different files on request !25
- Move read/write responsibility to model state !25
- States are now independent from models !25
- Output files changed its name scheme !25
- Configuration file changed its structure (see git diff on tests) !25
### Removed
- Coefficient mapper and its use on local operators !25
- Remove the operator map of objects within the model class !25
### Fixed
- Violation of the One-Definition-Rule due to external linkage of `LocalOperatorApply` lambdas !23
- Error on the gaussian equation used on comparisons !23

## [0.2.0] - 2020-02-20
### Added
- Code documentation
- Data Context concept for factories
- Factory concept for arbitrary object instantiation
- Add factories for finite element and finite element map
- Brief installation instructions
- Models can interpolate grid functions
- Grid utilities to recognize and mark tripes of entities
- A finite volume loca operator
- Grid function getters for external use
- Single domain executable
### Changed
- Move and rename header files
- Multidomain finite element was split into a multidomain and a dynamic power finite element
- Code DUNE dependencies are set to release 2.7 (See README.md)
- Other DUNE dependencies are set to COPASI/support/dune-copasi (See README.md)
- Executable is an optional build
- Library is optional and is split into smaller libraries
- Bump version utility updated to python3
### Fixed
- Dirichlet-Dirichlet condition at interfaces was being computed twice

## [0.1.0] - 2019-10-11
### Added
- [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
- Solver for reaction-diffusion systems in multiple compartments.

[Unreleased]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v0.2.0...master
[0.2.0]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v0.1.0...v0.2.0
[0.1.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/tags/v0.1.0
