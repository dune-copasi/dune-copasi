ARG BASE_IMAGE=debian:trixie

# ARG SETUP_BASE_IMAGE=registry.dune-project.org/docker/ci/${BASE_IMAGE}
ARG SETUP_BASE_IMAGE=${BASE_IMAGE}
ARG BUILD_BASE_IMAGE=setup-env
ARG PRODUCTION_BASE_IMAGE=${BASE_IMAGE}

# setup of dune dependencies
FROM ${SETUP_BASE_IMAGE} AS setup-env


RUN rm -f /etc/apt/apt.conf.d/docker-gzip-indexes \
  && rm -rf /var/lib/apt/lists/*
RUN export DEBIAN_FRONTEND=noninteractive; \
  apt-get update && apt-get dist-upgrade --no-install-recommends --yes \
  && apt-get install --no-install-recommends --yes \
  binutils \
  build-essential \
  ca-certificates \
  cmake \
  curl \
  dpkg \
  dpkg-dev \
  file \
  g++-13 \
  gcc-13 \
  gcovr \
  git \
  git-lfs \
  libfftw3-dev \
  libfftw3-mpi-dev \
  libgtest-dev \
  libmuparser-dev \
  libspdlog-dev \
  libsuitesparse-dev \
  libsuperlu-dev \
  libtbb-dev \
  libtiff-dev \
  mpi-default-bin \
  mpi-default-dev \
  ninja-build \
  pkg-config \
  && apt-get clean


ARG TOOLCHAIN=clang-6-17

ENV PATH=/duneci/install/bin:$PATH
ENV TERM=xterm-256color
ENV CMAKE_INSTALL_PREFIX=/duneci/install
ENV SETUP_DUNE_TESTTOOLS=ON
ENV DUNE_VENDOR_FMT=ON
# ENV DUNE_OPTS_FILE=/duneci/dune.opts
ENV DUNE_OPTS_FILE=/duneci/cmake-flags/dune-copasi.opts

COPY --chown=duneci ./dune-copasi.opts /duneci/cmake-flags/
COPY --chown=duneci ./.ci /duneci/modules/dune-copasi/.ci

RUN    ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain \
    && export PATH=/duneci/install/bin:$PATH
WORKDIR /duneci/modules
RUN mkdir -p /duneci/modules/dune-copasi/.ci
RUN ./dune-copasi/.ci/setup_dune $DUNE_OPTS_FILE

# build and install dune-copasi from the setup-env
FROM ${BUILD_BASE_IMAGE} AS build-env

ENV CPACK_GENERATORS=DEB
ENV CPACK_PACKAGE_DIRECTORY=/duneci/packages

WORKDIR /duneci/modules

# move source files into the image
COPY --chown=duneci ./ /duneci/modules/dune-copasi

# run installer
RUN ./dune-copasi/.ci/install $DUNE_OPTS_FILE

# tests installer
# RUN ./dune-copasi/.ci/test $DUNE_OPTS_FILE

# move results to a -lighter- production image
FROM ${PRODUCTION_BASE_IMAGE} AS production-env
LABEL maintainer="santiago.ospina@iwr.uni-heidelberg.de"

# get package from build-env and install it
COPY --from=build-env /duneci/packages/dune-copasi*Runtime.deb /packages/
WORKDIR /packages/
RUN rm -f /etc/apt/apt.conf.d/docker-gzip-indexes \
  && rm -rf /var/lib/apt/lists/*
RUN export DEBIAN_FRONTEND=noninteractive; \
  apt-get update && apt-get dist-upgrade --no-install-recommends --yes \
  && apt-get install --no-install-recommends --yes ./dune-copasi-*-Runtime.deb \
  && apt-get clean
RUN rm -rf /packages

# disable sudo user
RUN adduser --disabled-password --home /dunecopasi --uid 50000 dunecopasi
USER dunecopasi
WORKDIR /dunecopasi

# run version command and expect no error signal
RUN dune-copasi --version

# set default mout point to be /dunecopasi (same as workdir!)
VOLUME ["/dunecopasi"]
# run dune-copasi by default when running the image
ENTRYPOINT ["dune-copasi"]
