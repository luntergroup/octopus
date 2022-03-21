FROM ubuntu:impish

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Get dependencies
RUN apt-get -y update \
    && apt-get -y install \
        build-essential \
        libboost-all-dev \
        libgmp-dev \
        cmake \
        libhts-dev \
        python3-pip \
        git \
    && pip3 install distro

# Install Octopus
ARG threads=4
ARG architecture=haswell
COPY . /opt/octopus
RUN /opt/octopus/scripts/install.py \
    --threads $threads \
    --architecture $architecture

# Cleanup git - only needed during install for commit info
RUN apt-get purge -y git \
    && rm -r /opt/octopus/.git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/octopus/bin:${PATH}"

ENTRYPOINT ["octopus"]
