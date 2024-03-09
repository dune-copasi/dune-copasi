ARG BASE_IMAGE=debian:bookworm

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
  gcovr \
  git \
  git-lfs \
  libgtest-dev \
  ninja-build \
  pkg-config \
  python3 \
  python3-pip \
  python3-venv \
  nodejs \
  npm \
  && apt-get clean

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip3 install semver

ENV PATH=/duneci/install/bin:$PATH
ENV TERM=xterm-256color
ENV CMAKE_INSTALL_PREFIX=/duneci/install
ENV CMAKE_CXX_COMPILER=emcc
ENV CMAKE_C_COMPILER=emcc
ENV CMAKE_GENERATOR=Ninja
ENV CMAKE_TOOLCHAIN_FILE=/duneci/modules/emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake
ENV CMAKE_CROSSCOMPILING_EMULATOR=/duneci/modules/emsdk/node/16.20.0_64bit/bin/node
ENV DEFAULT_CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=em++ -DCMAKE_C_COMPILER=emcc -DBUILD_SHARED_LIBS=OFF -DCMAKE_CROSSCOMPILING_EMULATOR=${CMAKE_CROSSCOMPILING_EMULATOR} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE} -DCMAKE_CXX_FLAGS=-fexceptions"
ENV DUNE_OPTS_FILE=/duneci/cmake-flags/dune-copasi.opts

RUN useradd --no-log-init -u 50000 -g users --home /duneci duneci
WORKDIR /duneci/modules

RUN git clone https://github.com/emscripten-core/emsdk.git
RUN git clone https://github.com/oneapi-src/oneTBB.git
RUN git clone --depth 1 --branch v4.6.0 https://gitlab.com/libtiff/libtiff.git
RUN git clone --depth 1 --branch 9.1.0 https://github.com/fmtlib/fmt.git
RUN git clone --depth 1 --branch v1.11.0 https://github.com/gabime/spdlog.git

SHELL ["/bin/bash", "-c"]
# not working: 3.1.[52-54]
ENV EMSDK_VERSION=3.1.51
RUN ./emsdk/emsdk install ${EMSDK_VERSION}
RUN ./emsdk/emsdk activate ${EMSDK_VERSION}
ENV EMSDK_QUIET=1

RUN mkdir oneTBB/build \
    && source /duneci/modules/emsdk/emsdk_env.sh \
    && cmake oneTBB -B oneTBB/build -G Ninja $DEFAULT_CMAKE_FLAGS -DTBB_STRICT=OFF -DTBB_DISABLE_HWLOC_AUTOMATIC_SEARCH=ON -DTBB_EXAMPLES=OFF -DTBB_TEST=OFF
RUN cmake --build oneTBB/build
RUN cmake --install oneTBB/build
ENV TBB_ROOT=/duneci/install

COPY --chown=duneci ./dune-copasi.opts /duneci/cmake-flags/
COPY --chown=duneci ./.ci /duneci/modules/dune-copasi/.ci

ENV DUNE_PDELAB_ENABLE_TRACING=OFF
RUN    ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain \
    && export PATH=/duneci/install/bin:$PATH
RUN mkdir -p /duneci/modules/dune-copasi/.ci
RUN source /duneci/modules/emsdk/emsdk_env.sh \
    && ./dune-copasi/.ci/setup_dune $DUNE_OPTS_FILE

ENV dune-common_ROOT=/duneci/install
ENV dune-geometry_ROOT=/duneci/install
ENV dune-grid_ROOT=/duneci/install
ENV dune-istl_ROOT=/duneci/install
ENV dune-localfunctions_ROOT=/duneci/install
ENV dune-functions_ROOT=/duneci/install
ENV dune-uggrid_ROOT=/duneci/install
ENV dune-typetree_ROOT=/duneci/install
ENV dune-pdelab_ROOT=/duneci/install
ENV dune-multidomaingrid_ROOT=/duneci/install

RUN mkdir libtiff/build-tiff \
    && source /duneci/modules/emsdk/emsdk_env.sh \
    && cmake libtiff -B libtiff/build-tiff $DEFAULT_CMAKE_FLAGS -Dtiff-tests=OFF
RUN cmake --build libtiff/build-tiff
RUN cmake --install libtiff/build-tiff
ENV TIFF_ROOT=/duneci/install

RUN mkdir fmt/build \
    && source /duneci/modules/emsdk/emsdk_env.sh \
    && cmake fmt -B fmt/build $DEFAULT_CMAKE_FLAGS -DFMT_TEST=OFF
RUN cmake --build fmt/build
RUN cmake --install fmt/build
ENV fmt_ROOT=/duneci/install

RUN mkdir spdlog/build \
    && source /duneci/modules/emsdk/emsdk_env.sh \
    && cmake spdlog -B spdlog/build $DEFAULT_CMAKE_FLAGS -DSPDLOG_BUILD_TESTS=OFF -DSPDLOG_FMT_EXTERNAL=ON
RUN cmake --build spdlog/build
RUN cmake --install spdlog/build
ENV spdlg_ROOT=/duneci/install

RUN echo CMAKE_FLAGS+=\"\
 -Ddune-common_DIR=/duneci/install/lib/cmake/dune-common \
 -Ddune-geometry_DIR=/duneci/install/lib/cmake/dune-geometry \
 -Ddune-grid_DIR=/duneci/install/lib/cmake/dune-grid \
 -Ddune-istl_DIR=/duneci/install/lib/cmake/dune-istl \
 -Ddune-localfunctions_DIR=/duneci/install/lib/cmake/dune-localfunctions \
 -Ddune-functions_DIR=/duneci/install/lib/cmake/dune-functions \
 -Ddune-uggrid_DIR=/duneci/install/lib/cmake/dune-uggrid \
 -Ddune-typetree_DIR=/duneci/install/lib/cmake/dune-typetree \
 -Ddune-pdelab_DIR=/duneci/install/lib/cmake/dune-pdelab \
 -Ddune-multidomaingrid_DIR=/duneci/install/lib/cmake/dune-multidomaingrid \
 -Dspdlg_ROOT=/duneci/install \
 -Dfmt_ROOT=/duneci/install \
 -DTIFF_ROOT=/duneci/install \
 -DDUNE_COPASI_DISABLE_FETCH_PACKAGE_parafields=ON \
 -DDUNE_COPASI_DISABLE_FMT_STYLE=ON \
 -DCMAKE_EXE_LINKER_FLAGS=\'-sEXPORTED_RUNTIME_METHODS=ccall,cwrap,FS,callMain -sINVOKE_RUN=0 -sALLOW_MEMORY_GROWTH=1 -sEXPORT_ES6=1 -sMODULARIZE -sEXPORT_NAME=wasm -sENVIRONMENT=web -fexceptions\' \
 $DEFAULT_CMAKE_FLAGS \
 \" >> $DUNE_OPTS_FILE

# build and install dune-copasi from the setup-env
FROM ${BUILD_BASE_IMAGE} AS build-env

# move source files into the image
COPY --chown=duneci ./ /duneci/modules/dune-copasi

ENV CPACK_GENERATORS=TGZ
ENV CPACK_PACKAGE_DIRECTORY=/duneci/packages/

# run installer
RUN source /duneci/modules/emsdk/emsdk_env.sh \
    && ./dune-copasi/.ci/install $DUNE_OPTS_FILE
