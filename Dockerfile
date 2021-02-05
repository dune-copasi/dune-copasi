ARG BASE_IMAGE=debian:10
ARG DUNECI_PARALLEL=2

# setup of dune dependencies
FROM registry.dune-project.org/docker/ci/${BASE_IMAGE} AS setup-env

ARG TOOLCHAIN=clang-6-17

ENV DUNE_OPTIONS_FILE=/duneci/dune.opts
ENV PATH=/duneci/install/bin:$PATH
ENV TERM=xterm-256color
ENV CMAKE_INSTALL_PREFIX=/duneci/install

COPY --chown=duneci ./dune-copasi.opts /duneci/cmake-flags/dune-copasi.opts
RUN    ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain \
    && export PATH=/duneci/install/bin:$PATH
RUN    echo 'CMAKE_FLAGS+=" -DCMAKE_GENERATOR='"'"'Ninja'"'"' "' >> /duneci/cmake-flags/dune-copasi.opts \
    && echo 'CMAKE_FLAGS+=" -DDUNE_COPASI_SD_EXECUTABLE=ON"' >> /duneci/cmake-flags/dune-copasi.opts \

WORKDIR /duneci/modules
RUN mkdir -p /duneci/modules/dune-copasi/.ci
COPY --chown=duneci ./.ci /duneci/modules/dune-copasi/.ci
RUN RECURSIVE=ON SETUP_DUNE_TESTTOOLS=ON ./dune-copasi/.ci/setup_dune /duneci/dune.opts

# build and install dune-copasi from the setup-env
FROM setup-env AS build-env

ENV DUNE_OPTIONS_FILE=/duneci/dune.opts
ENV PATH=/duneci/install/bin:$PATH
ENV TERM=xterm-256color
ENV CMAKE_INSTALL_PREFIX=/duneci/install
ENV BUILD_TESING=ON

WORKDIR /duneci/modules
COPY --chown=duneci ./ /duneci/modules/dune-copasi
RUN ./dune-copasi/.ci/install /duneci/dune.opts

WORKDIR /duneci/modules/build-cmake/
RUN cpack -G DEB -B /duneci/packages CPackConfig.cmake

# move results to a -lighter- production  image
FROM ${BASE_IMAGE} AS production-env
LABEL maintainer="santiago.ospina@iwr.uni-heidelberg.de"

# get package from build-env and install it
COPY --from=build-env /duneci/packages/dune-copasi*Runtime.deb /packages/
WORKDIR /packages/
RUN rm -f /etc/apt/apt.conf.d/docker-gzip-indexes \
  && rm -rf /var/lib/apt/lists/*
RUN export DEBIAN_FRONTEND=noninteractive; \
  apt-get update && apt-get dist-upgrade --no-install-recommends --yes \
  && apt-get install --no-install-recommends --yes ./dune-copasi-*-Runtime.deb \
  && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN rm -rf /packages

# disable sudo user
RUN adduser --disabled-password --home /dunecopasi --uid 50000 dunecopasi
USER dunecopasi
WORKDIR /dunecopasi

# bring a test from build env
COPY --chown=dunecopasi --from=build-env /duneci/modules/dune-copasi/build-cmake/test/test_initial_triangles.ini /dunecopasi/
COPY --chown=dunecopasi --from=build-env /duneci/modules/dune-copasi/build-cmake/test/grids/square_triangles.msh /dunecopasi/grids/
# run and expect no signal if may runs
RUN dune-copasi-md test_initial_triangles.ini
RUN rm -rf *

# set default mout point to be /dunecopasi (same as workdir!)
VOLUME ["/dunecopasi"]
# run dune-copasi-md by default when running the image
ENTRYPOINT ["dune-copasi-md"]
