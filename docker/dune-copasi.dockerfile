ARG BASE_IMAGE=debian:bookworm

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
  doxygen \
  ca-certificates \
  cmake \
  codespell \
  curl \
  dpkg \
  dpkg-dev \
  file \
  gcovr \
  git \
  git-lfs \
  gnupg \
  libfftw3-dev \
  libfftw3-mpi-dev \
  libgtest-dev \
  libmuparser-dev \
  libscotchmetis-dev \
  libspdlog-dev \
  libsuitesparse-dev \
  libsuperlu-dev \
  libtbb-dev \
  libtiff-dev \
  mpi-default-bin \
  mpi-default-dev \
  ninja-build \
  pkg-config \
  python3 \
  python3-pip \
  python3-dev \
  python3-venv \
  python3-setuptools \
  software-properties-common \
  wget \
  zip \
  && apt-get clean

RUN  wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - \
  && echo "deb http://apt.llvm.org/bookworm/ llvm-toolchain-bookworm-18 main" > /etc/apt/sources.list.d/llvm.list \
  && export DEBIAN_FRONTEND=noninteractive; \
  apt-get update \
  && apt-get install --no-install-recommends --yes \
  clang-18 \
  clang-tidy-18 \
  && update-alternatives --install /usr/bin/clang clang /usr/bin/clang-18 100 \
  && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-18 100 \
  && update-alternatives --install /usr/bin/clang-tidy clang-tidy /usr/bin/clang-tidy-18 100 \
  && apt-get clean

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip3 install codechecker semver jschon

ENV PATH=/duneci/install/bin:$PATH
ENV TERM=xterm-256color
ENV CMAKE_INSTALL_PREFIX=/duneci/install
# ENV DUNE_OPTS_FILE=/duneci/dune.opts
ENV DUNE_OPTS_FILE=/duneci/cmake-flags/dune-copasi.opts

# Create user and group that the duneci will run under, and create base directory for builds
RUN adduser --disabled-password --home /duneci --group duneci duneci \
  && mkdir /builds
  && chown duneci:duneci /builds

USER duneci
COPY --chown=duneci ./dune-copasi.opts /duneci/cmake-flags/
COPY --chown=duneci ./.ci /duneci/modules/dune-copasi/.ci

ENV DUNE_ENABLE_PYTHONBINDINGS=OFF
ENV CMAKE_CXX_COMPILER=clang++
ENV CMAKE_C_COMPILER=clang
ENV CMAKE_GENERATOR="Ninja"
# RUN    ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain \
#     && export PATH=/duneci/install/bin:$PATH
WORKDIR /duneci/modules
RUN mkdir -p /duneci/modules/dune-copasi/.ci
RUN ./dune-copasi/.ci/setup_dune $DUNE_OPTS_FILE

# build and install dune-copasi from the setup-env
FROM ${BUILD_BASE_IMAGE} AS build-env

ENV CMAKE_CXX_COMPILER=clang++
ENV CMAKE_C_COMPILER=clang
ENV CMAKE_GENERATOR="Ninja"
ENV CPACK_GENERATORS=DEB
ENV CPACK_PACKAGE_DIRECTORY=/duneci/packages

WORKDIR /duneci/modules

# move source files into the image
COPY --chown=duneci ./ /duneci/modules/dune-copasi

# run installer
ENV DUNE_ENABLE_PYTHONBINDINGS=OFF
ENV DUNE_COPASI_DISABLE_FETCH_PACKAGE_parafields=OFF
RUN ./dune-copasi/.ci/install $DUNE_OPTS_FILE

# move results to a -lighter- production image
FROM ${PRODUCTION_BASE_IMAGE} AS production-env
LABEL maintainer="santiago@dune-project.org"

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

# test dune-copasi from the build-env
FROM build-env AS test-env
# tests installer
RUN ./dune-copasi/.ci/test $DUNE_OPTS_FILE

# run version command and expect no error signal
RUN dune-copasi --version

# run dune-copasi by default when running the image
ENTRYPOINT ["dune-copasi"]
