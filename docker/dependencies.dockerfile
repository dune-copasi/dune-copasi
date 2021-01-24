ARG BASE_IMAGE=registry.dune-project.org/docker/ci/ubuntu:18.04
FROM ${BASE_IMAGE}
ARG DUNECI_PARALLEL=2

ARG TOOLCHAIN=clang-6-17
ARG DUNE_OPTIONS_FILE=/duneci/dune.opts
ARG DUNECONTROL=/duneci/modules/dune-common/bin/dunecontrol

RUN ln -s /duneci/toolchains/${TOOLCHAIN} /duneci/toolchain

ENV PATH=/duneci/install/bin:$PATH

RUN echo 'CMAKE_FLAGS+=" -DDUNE_PYTHON_VIRTUALENV_SETUP=1 -DDUNE_PYTHON_VIRTUALENV_PATH=/duneci/install/dune-python-venv"' >> /duneci/cmake-flags/enable_virtualenv
RUN echo 'CMAKE_FLAGS+=" -DCMAKE_PREFIX_PATH:PATH=/duneci/install"' >> /duneci/cmake-flags/install_path
RUN echo 'CMAKE_FLAGS+=" -DCMAKE_INSTALL_PREFIX:PATH=/duneci/install"' >> /duneci/cmake-flags/install_path
COPY .ci/setup.sh /duneci/modules/setup.sh
RUN cd /duneci/modules/ && bash setup.sh
