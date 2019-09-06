ARG BASE_IMAGE=registry.dune-project.org/docker/ci/ubuntu:18.04
FROM ${BASE_IMAGE}
ARG DUNECI_PARALLEL=4

ARG TOOLCHAIN="clang-6-17"
RUN ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain

RUN duneci-install-module -b support/dune-copasi https://gitlab.dune-project.org/santiago.ospina/dune-common.git \
    && duneci-install-module --recursive https://gitlab.dune-project.org/staging/dune-logging.git \
    && duneci-install-module https://gitlab.dune-project.org/core/dune-geometry.git \
    && duneci-install-module https://gitlab.dune-project.org/staging/dune-uggrid.git \
    && duneci-install-module https://gitlab.dune-project.org/core/dune-grid.git \
    && duneci-install-module https://gitlab.dune-project.org/core/dune-istl.git \
    && duneci-install-module https://gitlab.dune-project.org/core/dune-localfunctions.git \
    && duneci-install-module -b support/dune-copasi https://gitlab.dune-project.org/santiago.ospina/dune-typetree.git \
    && duneci-install-module https://gitlab.dune-project.org/staging/dune-functions.git \
    && duneci-install-module -b support/dune-copasi https://gitlab.dune-project.org/santiago.ospina/dune-pdelab.git \
    && duneci-install-module -b support/dune-copasi https://gitlab.dune-project.org/santiago.ospina/dune-multidomaingrid.git \
    && duneci-install-module -b use-pdelab-dynamic-power-grid-function-space https://gitlab.dune-project.org/santiago.ospina/dune-copasi.git