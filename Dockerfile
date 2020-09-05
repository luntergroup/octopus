FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

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
COPY . /octopus
RUN octopus/scripts/install.py --dependencies --forests --threads 4

ENTRYPOINT ["octopus/bin/octopus"]
