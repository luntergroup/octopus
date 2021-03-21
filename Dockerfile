FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Get all apt dependencies
RUN apt-get -y update
RUN apt-get -y install software-properties-common
RUN apt-get install -y --no-install-recommends apt-utils
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get -y install \
    gcc g++ \
    build-essential \
    git \
    curl \
    python3-pip

RUN pip3 install distro

# Install Octopus
ARG threads=4
ARG architecture=haswell
COPY . /opt/octopus
RUN /opt/octopus/scripts/install.py \
    --dependencies \
    --forests \
    --threads $threads \
    --architecture $architecture

ENV PATH="/opt/octopus/bin:${PATH}"

ENTRYPOINT ["octopus"]
