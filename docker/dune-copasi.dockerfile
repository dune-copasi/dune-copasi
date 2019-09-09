ARG BASE_IMAGE
FROM ${BASE_IMAGE}

RUN duneci-install-module -b use-pdelab-dynamic-power-grid-function-space https://gitlab.dune-project.org/santiago.ospina/dune-copasi.git