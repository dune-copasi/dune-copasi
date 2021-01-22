ARG BASE_IMAGE
FROM ${BASE_IMAGE}
ARG DUNECI_PARALLEL=2
ARG BRANCH=master

RUN duneci-install-module -b ${BRANCH} https://gitlab.dune-project.org/copasi/dune-copasi.git
RUN dunecontrol --only=dune-copasi bexec cmake --build . -- install
RUN mkdir /duneci/dune-copasi && rm -rf /duneci/modules/

ENV PATH="${PATH}:/duneci/install/bin"
WORKDIR /duneci/dune-copasi
VOLUME ["/duneci/dune-copasi"]
