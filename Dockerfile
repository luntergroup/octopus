FROM ubuntu:impish

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Get all apt dependencies
RUN apt-get -y update
RUN apt-get -y install \
    build-essential \
    libboost-all-dev \
    libgmp-dev \
    cmake \
    libhts-dev \
    python3-pip

RUN pip3 install distro

# Install Octopus
ARG threads=4
ARG architecture=haswell
COPY . /opt/octopus
RUN /opt/octopus/scripts/install.py \
    --threads $threads \
    --architecture $architecture

ENV PATH="/opt/octopus/bin:${PATH}"

ENTRYPOINT ["octopus"]
