#!/bin/bash

echo "CC = ${CC}"
echo "CXX = ${CXX}"
echo "ASAN = ${ASAN}"
cmake --version;
${CC} --version;
${CXX} --version;

echo "COMMIT $(git rev-parse HEAD)"
echo "BRANCHES"
git branch
echo "----------"

if [ -n "${ASAN}" ]; then
  echo "Building with AddressSanitizer enabled..."
  SAN_FLAGS="-fsanitize=address -fsanitize=undefined -g -O1"
  export CXXFLAGS=${SAN_FLAGS}
  export LDFLAGS=${SAN_FLAGS}
fi

# Install apt-get dependencies
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    libicu-dev \
    curl \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcrypto++-dev \
    libboost-all-dev

# Install htslib
git clone https://github.com/samtools/htslib.git
cd htslib && autoheader && autoconf && ./configure && make && sudo make install
cd ..

./scripts/install.py --threads 1 --clean -c $CC -cxx $CXX
