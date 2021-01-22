ARG BUILD_IMAGE
ARG PRODUCTION_IMAGE

FROM ${BUILD_IMAGE} AS build-env

ARG BRANCH=master
ARG DUNECI_PARALLEL=2

RUN duneci-install-module -b ${BRANCH} https://gitlab.dune-project.org/copasi/dune-copasi.git
RUN dunecontrol --only=dune-copasi bexec cmake --build . -- install

FROM ${PRODUCTION_IMAGE} AS production-env
LABEL maintainer="santiago.ospina@iwr.uni-heidelberg.de"

COPY --from=0 /duneci/install /usr/local
RUN adduser --disabled-password --home /dunecopasi --uid 50000 dunecopasi

WORKDIR /dunecopasi
VOLUME ["/dunecopasi"]
USER dunecopasi
