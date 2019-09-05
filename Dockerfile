FROM ubuntu:latest

# Get all apt dependencies
RUN apt-get -y update
RUN apt-get -y install software-properties-common
RUN apt-get install -y --no-install-recommends apt-utils
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get -y install \
    gcc g++ \
    build-essential \
    make \
    wget \
    autotools-dev \
    libicu-dev \
    git \
    curl \
    libcurl4-openssl-dev \
    pkg-config \
    autoconf \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    openssl \
    libssl-dev \
    libcrypto++-dev

# Install Octopus
RUN git clone -b develop https://github.com/luntergroup/octopus.git
RUN octopus/scripts/install.py --install-dependencies --download-forests --threads 4

ENTRYPOINT ["octopus/bin/octopus"]
