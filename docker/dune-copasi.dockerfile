ARG BUILD_IMAGE
ARG PRODUCTION_IMAGE

FROM ${BUILD_IMAGE} AS build-env

ARG BRANCH=master
ARG DUNECI_PARALLEL=2

RUN    echo 'CMAKE_FLAGS=" -G '"'"'Unix Makefiles'"'"'"' >> /duneci/cmake-flags/production \
    && echo 'CMAKE_FLAGS+=" -DDUNE_USE_ONLY_STATIC_LIBS=ON"' >> /duneci/cmake-flags/production \
    && echo 'CMAKE_FLAGS+=" -DDUNE_COPASI_SD_EXECUTABLE=ON"' >> /duneci/cmake-flags/production \
    && echo 'CMAKE_FLAGS+=" -DDUNE_COPASI_MD_EXECUTABLE=ON"' >> /duneci/cmake-flags/production \
    && echo 'CMAKE_FLAGS+=" -DCMAKE_CXX_FLAGS='"'"'-fvisibility=hidden -fpic -static-libstdc++'"'"' "' /duneci/cmake-flags/production

RUN duneci-install-module -b ${BRANCH} https://gitlab.dune-project.org/copasi/dune-copasi.git
RUN dunecontrol --opts=/duneci/dune.opts --only=dune-copasi bexec cmake --build . -- install

FROM ${PRODUCTION_IMAGE} AS production-env
LABEL maintainer="santiago.ospina@iwr.uni-heidelberg.de"

COPY --from=0 /duneci/install /usr/local
RUN adduser --disabled-password --home /dunecopasi --uid 50000 dunecopasi

WORKDIR /dunecopasi
VOLUME ["/dunecopasi"]
USER dunecopasi
