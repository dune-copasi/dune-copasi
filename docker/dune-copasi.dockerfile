ARG BASE_IMAGE
FROM ${BASE_IMAGE}
ARG DUNECI_PARALLEL=2

RUN duneci-install-module https://gitlab.dune-project.org/santiago.ospina/dune-copasi.git