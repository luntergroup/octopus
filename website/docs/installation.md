---
id: installation
title: Installation
---

Octopus can be built and installed on most Unix based systems (e.g. Linux and MacOS). Windows has not been tested.

The recommend way to install Octopus for most users is:

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git
$ octopus/scripts/install.py --dependencies --forests
```

You can then optionally add `octopus` to your `PATH`:

```shell
$ echo 'export PATH='$(pwd)'/octopus/bin:$PATH' >> ~/.bash_profile
$ source ~/.bash_profile
```

Then check the installation was successful:

```shell
$ octopus --version
```

## Requirements

* A [C++14](https://isocpp.org/wiki/faq/cpp14) compiler and compatibility standard library. Either [GCC](https://gcc.gnu.org) (version >= 9.3) or [Clang](https://clang.llvm.org) (version >= 11.0) are recommended.
* [Git](https://git-scm.com) version >= 2.5
* [Boost](https://www.boost.org) version >= 1.65
* [htslib](https://github.com/samtools/htslib) version >= 1.4; version != 1.12
* [GMP](https://gmplib.org) version >= 5.1.0
* [CMake](https://cmake.org) version >= 3.9
* Optional:
    * [Python](https://www.python.org) version >= 3 plus the [distro](https://pypi.org/project/distro/) package

:::important

Octopus uses [SIMD](https://en.wikipedia.org/wiki/SIMD) instructions for performance reasons. The instruction set used (minimum [SSE2](https://en.wikipedia.org/wiki/SSE2)) is built statically, so if you compile with [AVX2](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2), you won't be able to use the resulting binary on machines that doesn't support AVX2.

:::

## Python

First clone the git repository in your preferred directory:

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git && cd octopus
```

The easiest way to install Octopus from source is with the Python3 installer script. To see the options available to this script run `scripts/install.py --help`.

If all the requirements are accessible on your `PATH` then simply run

```shell
$ scripts/install.py
```

otherwise you can specify paths to each dependency, for example, to set the compiler you'd use 

```shell
$ scripts/install.py --cxx_compiler /path/to/cpp/compiler
```

By default, this installs to `/bin` relative to where octopus is installed. To install to a different location (e.g. `/usr/local/bin`) use:

```shell
$ scripts/install.py --prefix /user/local/bin
```

You can also request all dependencies to be installed locally:

```shell
$ scripts/install.py --dependencies
```

:::tip

If a build isn't working after an update then try adding `--clean` to the install command.

:::

### Setting the build architecture

By default, the binary is optimised for the build machine architecture. If you need to run Octopus on another machine with a different architecture then use the `--architecture` option:

```shell
$ scripts/install.py --architecture haswell
```

This is passed to the [-march](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html) compiler option. 

## CMake

If Python3 isn't available, Octopus can be installed directly with [CMake](https://cmake.org):

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git
$ cd octopus/build
$ cmake .. && make install
```

CMake will try to find a suitable compiler on your system, if you'd like you use a specific compiler use the `-D` option, for example:

```shell
$ cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ ..
```

## Docker

Pre-built Docker images are available on [DockerHub](https://hub.docker.com/r/dancooke/octopus):

```shell
$ docker pull dancooke/octopus
$ docker run dancooke/octopus -h
```

:::important

The Octopus images on DockerHub are built to a [Haswell](https://en.wikipedia.org/wiki/Haswell_(microarchitecture)) architecture. This means that they will only work on  Haswell (with AVX2) or newer machines.

:::

You can also build a new image from the Dockerfile:

```shell
$ git clone https://github.com/luntergroup/octopus.git && cd octopus
$ docker build -t octopus .
$ docker run octopus -h
```

This is especially useful if you need to build to a specific architecture:

```shell
$ docker build -t octopus --build-args architecture=sandybridge .
```

## Singularity

To build a [Singularity](https://singularity.hpcng.org) container directly from the DockerHub images use

```shell
$ singularity build octopus.sif docker://dancooke/octopus
```

## Conda

Octopus is available [pre-built for Linux](https://anaconda.org/bioconda/octopus) as part of [Bioconda](https://bioconda.github.io/):

```shell
$ conda install -c bioconda octopus
```

:::important

The Octopus package on Bioconda is built to a [Haswell](https://en.wikipedia.org/wiki/Haswell_(microarchitecture)) architecture. This means that it will only work on  Haswell (with AVX2) or newer machines. If you need another architecture then consider using [conda-build](https://docs.conda.io/projects/conda-build/en/latest/).

:::