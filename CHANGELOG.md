# Changelog

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

## [Unreleased] ([git-diff][Unreleased-diff])

## [2.1.0] ([git-diff][2.1.0-diff]) - 2024-11-25

### Added

 - Add experimental support for moving domains !205
 - Enable super builds for simpler installation instructions !230
 - Numerical jacobians !204

### Changed

 - Update ExprTk to 0.0.3 !237
 - Improve documentation !210 !211 !212 !213 !214 !215 !216 !217 !218 !219 !220 !221 !222 !223 !224 !225 !226 !227 !228

### Fixed

 - Wrong configuration of higher order polynomial degrees !236
 - Memory leak on the UMFPack solver !235
 - Outflow derivatives set its prefix incorrectly !229
 - Improve muParser search with CMake !203 !232

## [2.0.1] ([git-diff][2.0.1-diff]) - 2024-04-22

### Fixed

 - Fix memory errors on Blocked Jacobi !197
 - Remove Muparser invalid assertions !196

## [2.0.0] ([git-diff][2.0.0-diff]) - 2024-03-26

_**Note**: This version is a complete re-write of the library based on a custom version of PDELab tailored for this module. Thus, this changelog does not make much sense in comparison to version [1.1.0]. Instead, the following entries are written with respect to the first refactor from [1.1.0] made on !83._

### Refactor

 - Allows different parser backends: [muparser](https://beltoforion.de/en/muparser/), [exprtk](https://github.com/ArashPartow/exprtk), and [SymEngine](https://github.com/symengine/symengine) !83
 - Adds constants, (lambda) function, and [random fields](https://github.com/parafields/parafields-core) definitions within the parsers !83
 - Re-implementation of local operator: huge performance improvements !83
 - Adds tensor cross-diffusion, storage, velocity, and external boundary terms !83
 - Allows for volume coupling between compartments !83
 - Splits the subdomain and the compartment concepts. Now compartments with no entities or no components are perfectly possible !83
 - Uses new basis functions from dune-pdelab that allows native usage of multi-domains (branch: features/dune-assembler/main) !83
 - Extends compatibility of tiff images to different bit sizes !83
 - Adds monitoring (info, warning and error) of variables with new generic reduce operators !83
 - Switches logging from dune-logging to [spdlog](https://github.com/gabime/spdlog) !83
 - Designed to be thread-safe, although is not entirely yet implemented !83
 - Adds tracing capabilities with [perfetto](https://perfetto.dev/) !83
 - Allows Selection of between direct and sparse solvers, and matrix free operators !83
 - Allows to print matrix layout in a SVG file !83
 - Unifies the executable for different dimensions and different degrees of freedom layouts !83
 - Improves command line interface !83

### Added
 - Docusaurus now upgrades latest `canary` wasm binary automatically !186
 - Hierarchical ISTL solvers and preconditioners in a dynamic registry !156
 - Initial web interface for Wasm executable !170
 - Deploy NPM packages to GitLab registry !179
 - Tutorial and test on Cardiac Electrophysiology simulations !157
 - Grid cell data !175
 - Push CI artifacts to the GitLab registry !169 !173
 - Test more complicated cases !158
 - Job to build Wasm binaries in the CI !150
 - Math expressions with only numbers as scientific notation don't need a parser !154
 - Allow underscores in parser function arguments !152
 - Make multi-threading optional !149
 - Compile 2D and 3D in GitHub Actions !141
 - Customization point on grid axis names !139
 - Diagnostic information on SymEngine compilation failure !134
 - Test on pattern generation !132
 - Multi-threaded assembly !98
 - Allow arbitrary number of functions in ExprTk !148
 - Codequality check on the CI !123
 - LLVM visitor for SymEngine !121
 - Throw meaningful error message when config file is empty !119
 - Diagnostic information on ExprTk compilation failure !115
 - Command line help information !107 !120 !138
 - Interpolation for 1D functions !102
 - Initial value to reuction operations !100
 - Code spelling to the CI !88
 - Helper constraints for function spaces !87
 - Support libc++ and reduce minimum required standard from C++23 to C++20 !84
 - Show error when config file does not exist !86
 - ~~Doxygen documentation is now deployed to Netlify !73~~ (reverted on !176)
### Changed
 - Update docusaurus documentation to latest changes !176 !178 !88
 - Math parsers now can only take move only functions !187
 - Reduce operations now receive one less argument !184
 - Reorganize files and add `DiffusionReaction` namespace !183
 - Use grid leaf view instead of level 0 grid view for grid cell data !177
 - Author email address !181
 - Improve support for Wasm binaries !170 !182
 - Move configuration options to a JSON file and use schema validator on it !174
 - Use Debian Bookworm in the CI !171
 - Conditionally test possion config if ExprTk is available !168
 - Use (faster) integer based multi-domain grid !165
 - Simplify CMake usage !161 !162 !163
 - Improve version handler !159
 - Improve local basis cache to also handle intersections !145
 - Improve log error on failure !136
 - How to use constraints !135
 - Update docusaurus to version `3.0.0` !124
 - Cleanup of the build system files !118
 - Cleanup the CI installation scripts !110
 - Make options consistent with respect to file inputs `--*.path=/path/to/file` !109
 - Use `PDELab::Execution` to express (possible) concurrency !108
 - Use newer commit on `parafields-core` !105
 - Improve implementation of muParser and ExprTk parsers !96
 - Update docusaurus to version `2.0.0-alpha.75` !72
### Fixed
 - Missing CMake installation of Hierarchical ISTL solvers (from !156) !183
 - Boundary constraints were handled incorrectly !167
 - Coloring option for multi-threading was swapped with micro-locks !166
 - Remove wrong definition of maximum number of ExpTk functions !130
 - Wrong definition for SuiteSparse library !144
 - Missing user in Dockerfile !153 !155
 - ~~Freeze `dune-grid` dependency !147~~ (reverted by !164)
 - Use the keyword `no_value` to identify special values in the parsers !140
 - Inheritance of `parser_context.parser_type` !133
 - Entries per row estimate for pattern creation !127
 - Use `NDEBUG` compiler definition on release build !131
 - Compatibility issue with {fmt} 10.2 !125
 - Duplicated output on reduction algorithms !117
 - Testing in the CI !114
 - Correct construction on parser contexts !116
 - Compatibility with XCode and Msys !111
 - Consistent use of `struct` across different headers !104
 - Unused function keywords on SymEngine where throwing unnecessary errors !101
 - Deprecation warning on SVG writer !97
 - Compile final executable with different parsers !95
 - CMake problems on target installation !94
 - Support for {fmt} >= 9.0.0 !90 !91 !93 !129
 - Docker image now uses `dune-copasi` executable instead of `dune-copasi-[sd|md]` !91
 - CI jobs now passes since the re-write !89
 - Solver on linear problems was wrongly reused !85
### Removed
 - Unused header files from version [1.1.0] !183
 - Dependency on `dune-testtools` !122
 - Inherited `parafield-core` tests !112
 - Fallback for `std::filesystem` !94
###  Security
 - Update `follow-redirects` in web documentation !126

## [1.1.0] ([git-diff][1.1.0-diff]) - 2021-04-29
### Added
 - Dune options file receives `CMAKE_OPTIONS` and `MAKE_OPTIONS` !60
 - Parameter to use max norm on newton steps !65
 - Parameter to accept best solution in linear search on newton steps !65
 - Event timestepper !64
 - Preprocessor macros containing version information !67
### Changed
 - TIFF helper is compiled in the dunecopasi library !62
 - Use `NewtonMethod` instead of `Newton` class from PDELab !65
### Fixed
 - Executables can be compiled without the library !60
 - All dune cmake flags are now transitively passed to dune-copasi targets !60
 - System tests on main executables were ignored !66
 - Snap to final time now avoids very small timesteps !61
 - CMake config file was installed on the wrong path !68

## [1.0.0] ([git-diff][1.0.0-diff]) - 2021-02-11
### Added
- Custom membrane flux !30
- Skip intersection methods in local operators !36
- Parameters to control Newton's method !30
- Stepper interface to snap solution to a specific time !45
- Provide help message when executables are not used correctly !49
- Produce Debian Packages on the CI !49
- Installation is now divided on three components: `Runtime|Library|Development` !49
- Define a recommended dune options file [`dune-copasi.opts`](dune-copasi.opts) !49
- Versioned documentation !49
- Dependency on `pkg-config` !52
### Changed
- TIFF images are clamped instead of zeroed when evaluated outside its domain !47
- Produce versioned documentation !49
- Travis and Appveyor are replaced for GitHub Actions on `Linux|MacOS|Windows` !49
- Increment CMake required version to 3.13 !49
- Installed CMake project is now consumable by other CMake projects !49
- Vendor GHC with CMake instead of git submodules !49
- Increase vendored GHC version to 10.5 !49 & !51
- Force GHC usage whenever c++ filesystem cannot be found, else, optional !49
- CI scripts are improved to simplify usage !49
- Simplify build and usage on docker containers !49
- Drastically reduce size of final docker container !49
- Improve performance on cases with no interaction between all species !43
- `dune-logging` and `dune-multidomaingrid` no longer require a COPASI namespace fork !53
### Removed
- Automatic flux between compartment components with same name !30
- Jacobian operator is managed by the stepper instead of the model !39
- Jacobian tests are improved and extended to cover more cases !39
- CMake module to find muparser, use `pkg-config` instead !52
### Fixed
- Finished and documented Installation procedure !49
- Wrong rounding on y pixels on TIFF file reads !47
- Improper vector allocation in multidomain intersections !39
- Final timestep reduction to reach `end_time` failed in adaptive stepper !43
- Infinite loop when final adaptive step failed !45
- Binary executables are now installed !49
- Support for shared libraries !54

## [0.3.0] ([git-diff][0.3.0-diff]) - 2020-10-07
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

## [0.2.0] ([git-diff][0.2.0-diff]) - 2020-02-20
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

[Unreleased-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v2.0.1...master
[2.0.1-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v2.0.1...v2.0.1
[2.0.0-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v1.1.0...v2.0.0
[1.1.0-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v1.0.0...v1.1.0
[1.0.0-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v0.3.0...v1.0.0
[0.3.0-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v0.2.0...v0.3.0
[0.2.0-diff]: https://gitlab.dune-project.org/copasi/dune-copasi/compare/v0.1.0...v0.2.0

[Unreleased]: https://gitlab.dune-project.org/copasi/dune-copasi/-/tree/master

[2.0.1]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v2.0.1
[2.0.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v2.0.0
[1.1.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v1.1.0
[1.0.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v1.0.0
[0.3.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v0.3.0
[0.2.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v0.2.0
[0.1.0]: https://gitlab.dune-project.org/copasi/dune-copasi/-/releases/v0.1.0
