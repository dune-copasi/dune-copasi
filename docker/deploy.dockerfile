ARG PRODUCTION_BASE_IMAGE=${BASE_IMAGE}

# move results to a -lighter- production image
FROM ${PRODUCTION_BASE_IMAGE} AS production-env
LABEL maintainer="santiago.ospina@iwr.uni-heidelberg.de"

# get package from host 'packages' directory and install it
COPY ./packages /packages/
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

# run help and expect no error signal
RUN dune-copasi-md --help
RUN dune-copasi-sd --help

# set default mout point to be /dunecopasi (same as workdir!)
VOLUME ["/dunecopasi"]
# run dune-copasi-md by default when running the image
ENTRYPOINT ["dune-copasi-md"]
