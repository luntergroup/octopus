FROM ubuntu:17.10

# Get all apt dependencies
RUN apt-get -y update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get -y install \
    gcc \
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

# Install GCC 7
RUN apt-get -y install software-properties-common python-software-properties
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get -y install \
    gcc-7 \
    libstdc++6

# Install CMake
WORKDIR /tmp
RUN wget https://cmake.org/files/v3.10/cmake-3.10.2.tar.gz
RUN tar -xzvf cmake-3.10.2.tar.gz
WORKDIR /tmp/cmake-3.10.2
RUN ./bootstrap --prefix=/usr/local
RUN make -j2
RUN make install

# Install Boost
WORKDIR /tmp
RUN wget -O boost_1_66_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.66.0/boost_1_66_0.tar.gz/download
RUN tar xzvf boost_1_66_0.tar.gz
WORKDIR /tmp/boost_1_66_0
RUN ./bootstrap.sh --prefix=/usr/local --without-libraries=python,mpi
RUN ./b2 -j2 cxxflags="-std=c++11"
RUN ./b2 install

# Install htslib
WORKDIR /tmp
RUN git clone -b master https://github.com/samtools/htslib.git
WORKDIR /tmp/htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make -j2
RUN make install

# Install Octopus
WORKDIR /tmp
RUN git clone -b master https://github.com/luntergroup/octopus.git
WORKDIR /tmp/octopus
RUN ./install.py --root --threads=2

WORKDIR /home
