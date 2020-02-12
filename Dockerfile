# Build this image using the the file called build.sh
ARG BASE_IMAGE=ubuntu:18.04
FROM ${BASE_IMAGE}

ARG BUILD_TYPE=Release
ARG OMP_NUM_THREADS=16

# The following environment variables help with the Python 2 vs 3 option while
# we still support both of them.  In order for this approach to work, we need
# the RUN command to use bash instead of sh for variable indirection.
SHELL ["/bin/bash", "-c"]

# -----------------------------------------------------------------------------
# Timezone handling for tzdata
# -----------------------------------------------------------------------------

ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# -----------------------------------------------------------------------------
# System development package
# -----------------------------------------------------------------------------

RUN apt-get update && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        ca-certificates \
        chrpath \
        curl \
        gfortran \
        git \
        libtool \
        openssl \
        openmpi-bin \
        libopenmpi-dev \
        tcl-dev tk-dev libopenblas-dev liblapack-dev openssh-server \
        pkg-config && \
        rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# Non privilege user: ubuntu
# -----------------------------------------------------------------------------

# Create a non-root user
RUN groupadd ubuntu && \
    useradd -g ubuntu -d /home/ubuntu ubuntu && \
    mkdir /home/ubuntu && chown -R ubuntu:ubuntu /home/ubuntu
USER ubuntu

# -----------------------------------------------------------------------------
# EcoSlim
# -----------------------------------------------------------------------------

RUN cd && \
    git clone --recursive https://github.com/reedmaxwell/EcoSLIM.git

WORKDIR /home/ubuntu/EcoSLIM

# Build
RUN /usr/bin/make

ENV PATH /home/ubuntu/EcoSLIM:$PATH

USER ubuntu
ENTRYPOINT ["/bin/bash"]

