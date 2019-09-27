ARG BASE_IMAGE
FROM ${BASE_IMAGE}
ARG DUNECI_PARALLEL=2
ARG BRANCH=master

RUN duneci-install-module -b ${BRANCH} https://gitlab.dune-project.org/copasi/dune-copasi.git