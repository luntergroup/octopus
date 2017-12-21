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
RUN wget https://cmake.org/files/v3.9/cmake-3.9.6.tar.gz
RUN tar -xzvf cmake-3.9.6.tar.gz
WORKDIR /tmp/cmake-3.9.6
RUN ./bootstrap --prefix=/usr/local
RUN make -j4
RUN make install

# Install Boost
WORKDIR /tmp
RUN wget -O boost_1_65_1.tar.gz http://sourceforge.net/projects/boost/files/boost/1.65.1/boost_1_65_1.tar.gz/download
RUN tar xzvf boost_1_65_1.tar.gz
WORKDIR /tmp/boost_1_65_1
RUN ./bootstrap.sh --prefix=/usr/local --without-libraries=python,mpi
RUN ./b2 -j4 cxxflags="-std=c++11"
RUN ./b2 install

# Install htslib
WORKDIR /tmp
RUN git clone https://github.com/samtools/htslib.git
WORKDIR /tmp/htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make -j4
RUN make install

# Install Octopus
WORKDIR /tmp
RUN git clone https://github.com/luntergroup/octopus.git
WORKDIR /tmp/octopus
RUN ./install.py --root --threads=2

WORKDIR /home
