FROM ubuntu:latest

# Get all apt dependencies
RUN apt-get -y update
RUN apt-get -y install software-properties-common
RUN apt-get install -y --no-install-recommends apt-utils
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
RUN apt-get -y install \
    gcc g++ \
    build-essential \
    git \
    curl 

# Install Octopus
COPY . /octopus
RUN octopus/scripts/install.py --install-dependencies --download-forests --threads 4

ENTRYPOINT ["octopus/bin/octopus"]
