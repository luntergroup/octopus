FROM ubuntu:latest

# Get all apt dependencies
RUN apt-get -y update
RUN apt-get -y install software-properties-common
RUN apt-get install -y --no-install-recommends apt-utils
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get -y install \
    gcc-8 g++-8 \
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
ENV CC gcc-8
ENV CXX g++-8

# Install CMake
WORKDIR /tmp
RUN wget https://cmake.org/files/v3.11/cmake-3.11.4.tar.gz
RUN tar -xzvf cmake-3.11.4.tar.gz
WORKDIR /tmp/cmake-3.11.4
RUN ./bootstrap --prefix=/usr/local
RUN make -j2
RUN make install

# Install Boost
WORKDIR /tmp
RUN wget -O boost_1_68_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.68.0/boost_1_68_0.tar.gz/download
RUN tar xzvf boost_1_68_0.tar.gz
WORKDIR /tmp/boost_1_68_0
RUN ./bootstrap.sh --prefix=/usr/local --without-libraries=python,mpi
RUN ./b2 -j2 toolset=gcc-8 cxxflags="-std=c++14"
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
RUN git clone -b develop https://github.com/luntergroup/octopus.git
WORKDIR /tmp/octopus
RUN ./scripts/install.py --root --threads=2
RUN ldconfig

WORKDIR /home