ARG BASE_IMAGE=debian:10
ARG DUNECI_PARALLEL=2

# setup of dune dependencies
FROM registry.dune-project.org/docker/ci/${BASE_IMAGE} AS setup-env

ARG TOOLCHAIN=clang-6-17

ENV DUNE_OPTIONS_FILE=/duneci/dune.opts
ENV PATH=/duneci/install/bin:$PATH

RUN    ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain \
    && export PATH=/duneci/install/bin:$PATH
RUN    echo 'CMAKE_FLAGS+=" -DDUNE_PYTHON_VIRTUALENV_SETUP=1"' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DDUNE_PYTHON_VIRTUALENV_PATH=/duneci/install/dune-python-venv"' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DCMAKE_PREFIX_PATH:PATH=/duneci/install"' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DCMAKE_INSTALL_PREFIX:PATH=/duneci/install"' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DCMAKE_GENERATOR="Ninja' >> /duneci/cmake-flags/dune-copasi
WORKDIR /duneci/modules
RUN mkdir -p /duneci/modules/dune-copasi/.ci
COPY --chown=duneci ./.ci /duneci/modules/dune-copasi/.ci
RUN bash dune-copasi/.ci/setup.sh

# build and install dune-copasi from the setup-env
FROM setup-env AS build-env

ENV DUNE_OPTIONS_FILE=/duneci/dune.opts
ENV PATH=/duneci/install/bin:$PATH

RUN    echo 'CMAKE_FLAGS+=" -DDUNE_COPASI_SD_EXECUTABLE=ON"' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DDUNE_COPASI_MD_EXECUTABLE=ON"' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DCMAKE_CXX_FLAGS_RELEASE='"'"'-O3 -fvisibility=hidden -fpic -static-libstdc++'"'"' "' >> /duneci/cmake-flags/dune-copasi \
    && echo 'CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"' >> /duneci/cmake-flags/dune-copasi

WORKDIR /duneci/modules
COPY --chown=duneci ./ /duneci/modules/dune-copasi
RUN bash dune-copasi/.ci/build.sh

# move results to a -lighter- production  image
FROM ${BASE_IMAGE} AS production-env
LABEL maintainer="santiago.ospina@iwr.uni-heidelberg.de"

COPY --from=build-env /duneci/install /usr/local

RUN rm -f /etc/apt/apt.conf.d/docker-gzip-indexes \
  && rm -rf /var/lib/apt/lists/*
RUN export DEBIAN_FRONTEND=noninteractive; \
  apt-get update && apt-get dist-upgrade --no-install-recommends --yes \
  && apt-get install --no-install-recommends --yes \
  libscotchparmetis-dev \
  libldl2 \
  libspqr2 \
  libumfpack5 \
  libarpack++2c2a \
  libsuperlu5 \
  libgmpxx4ldbl \
  libopenblas-base \
  libtiff5 \
  libmuparser2v5 \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN adduser --disabled-password --home /dunecopasi --uid 50000 dunecopasi

WORKDIR /dunecopasi
VOLUME ["/dunecopasi"]
USER dunecopasi
