[![Build Status](https://gitlab.dune-project.org/copasi/dune-copasi/badges/master/pipeline.svg)](https://gitlab.dune-project.org/copasi/dune-copasi/pipelines)
[![Build Status](https://github.com/dune-copasi/dune-copasi/workflows/CI%20Builds/badge.svg?branch=master)](https://github.com/dune-copasi/dune-copasi/actions?query=branch%3Amaster+)

## How does the CI works?

We have 2 Continous Integration services:
  - [GitLab](https://docs.gitlab.com/ee/ci/)
  - [GitHub](https://github.com/dune-copasi/dune-copasi/actions)

For this, we have set up a main repository and a mirror:

  - Main repository: https://gitlab.dune-project.org/copasi/dune-copasi
  - Mirror repository: https://github.com/dune-copasi/dune-copasi

The main idea here is that they, the different CI, follow almost the same
instructions in the different stages. In all the cases, the scripts expect
configuration options defined in the main directory `dune-copasi.opts`.

Their files are defined in the default places:

  - [../.gitlab-ci.yml](../.gitlab-ci.yml)
  - [../.github/workflow/ci.yml](../.github/workflow/ci.yml)

### Stage: Setup Precopiled Dependecies

In the case of GitHub, we use the precopiled libraries provided by the
[`sme_deps_common`](https://github.com/spatial-model-editor/sme_deps_common).
The script [`setup_static_deps`](setup_static_deps) obtains the target OS
from the `dune-copasi.opts`, downloads the respective files, and expand them
on the cmake install prefix.

```bash
./setup_static_deps $PDW/../dune-copasi.opts
```

### Stage: Setup Dune Dependecies

In this stage, we configure and install the dune dependencies that `dune-copasi`
requires. In both cases, the CIs use the script [`setup_dune`](setup_dune), but
in order to generate faster results in GitLab, this stage is only updated for the
`master` branch by the docker file [`Dockerfile`](../docker/dune-copasi.dockerfile)
with target `setup-env`.

```bash
./setup_dune $PDW/../dune-copasi.opts
```

### Stage: Build, Install, and Package

This stage builds, installs and optionally packages `dune-copasi` following the
[`install`](install) script. Similar to the others, this script expects the options
file which in turs run `cmake` and its targets. If the environmental variables
`CPACK_GENERATORS` and `CPACK_PACKAGE_DIRECTORY` are set, the script will package the
library with the given generators (e.g. debian packages).

```bash
./install $PDW/../dune-copasi.opts
```

### Stage: Testing

This stage builds and runs the unit and system tests by using the
[`test`](test) script. `ctests` runs all tests with the tag
`unit` and `system`.


```bash
./test $PDW/../dune-copasi.opts
```

## CI Link & Badge
For getting them replace `<branch>` for the specific branch you want:

### GitLab CI
  - Link: https://gitlab.dune-project.org/copasi/dune-copasi/pipelines
  - Badge: https://gitlab.dune-project.org/copasi/dune-copasi/badges/`<branch>`/pipeline.svg
### GitHub Actions
  - Link: https://github.com/dune-copasi/dune-copasi/actions?query=branch%3A`<branch>`+
  - Badge: https://github.com/dune-copasi/dune-copasi/workflows/CI%20Builds/badge.svg?branch=`<branch>`
