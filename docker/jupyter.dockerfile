# This Dockerfile tries to be compatible with
# https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile
FROM jupyter/base-notebook:584f43f06586

USER root
RUN apt update && \
    apt install --no-install-recommends --yes \
      libgl1-mesa-glx && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*
USER ${NB_USER}

# install conda and its dependencies
RUN conda install -c conda-forge \
        cmake \
        git \
        cxxopts \
        fortran-compiler \
        gcc_linux-64 \
        gnuplot \
        make \
        suitesparse \
        superlu \
        muparser \
        libtiff \
        pkg-config \
        xeus-cling && \
    conda clean -a -q -y

WORKDIR ${HOME}/dune-modules

# set variables for the options file
ENV CMAKE_C_COMPILER="'/opt/conda/bin/x86_64-conda-linux-gnu-gcc'"
ENV CMAKE_CXX_COMPILER="'/opt/conda/bin/x86_64-conda-linux-gnu-g++'"
ENV CMAKE_CXX_FLAGS="''"
ENV CMAKE_OPTIONS='-DBUILD_SHARED_LIBS=ON'
ENV CMAKE_INSTALL_PREFIX="/opt/conda"

# copy options file
COPY --chown=${NB_UID} ./dune-copasi.opts ${HOME}/dune-modules/dune-copasi/dune-copasi.opts

# copy installation scripts
COPY --chown=${NB_UID} ./.ci/ ${HOME}/dune-modules/dune-copasi/.ci
RUN ./dune-copasi/.ci/setup_dune ${HOME}/dune-modules/dune-copasi/dune-copasi.opts

# copy main repository
COPY --chown=${NB_UID} ./ ${HOME}/dune-modules/dune-copasi

# set up variables for dune-copasi configuration
ENV DUNE_COPASI_SD_EXECUTABLE='OFF'
ENV DUNE_COPASI_MD_EXECUTABLE='OFF'
ENV VERBOSE=1

# configure and build dune-copasi
RUN ./dune-copasi/dune-copasi.opts
RUN ./dune-copasi/.ci/install ${HOME}/dune-modules/dune-copasi/dune-copasi.opts

# configure jupyter notebooks
RUN mkdir dune-copasi/notebooks-build \
    && cd dune-copasi/notebooks-build \
    && cmake -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} \
             -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} \
             -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} \
             -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX} \
             ../notebooks \
    && make \
    && make install

# make notebooks the entry point
WORKDIR ${HOME}/dune-modules/dune-copasi/notebooks

# We enable JupyerLab!
ENV JUPYTER_ENABLE_LAB=yes
