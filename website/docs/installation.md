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

## Requirements



## Python

First clone the git repository in your preferred directory:

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git && cd octopus
```

The easiest way to install octopus from source is with the Python3 install script. If your default compiler satisfies the minimum requirements just execute:

```shell
$ ./scripts/install.py
```

otherwise explicitly specify the compiler to use:

```shell
$ ./scripts/install.py --cxx_compiler /path/to/cpp/compiler # or just the compiler name if on your PATH
```

For example, if the requirement instructions above were used:

```shell
$ ./scripts/install.py --cxx_compiler clang++-4.0
```

On some systems, you may also need to specify a C compiler which is the same version as your C++ compiler, otherwise you'll get link errors. This can be done with the `--c_compiler` option:

```shell
$ ./scripts/install.py -cxx g++-7 -c gcc-7 
```

By default this installs to `/bin` relative to where you installed octopus. To install to a different location (e.g. `/usr/local/bin`) use:

```shell
$ ./scripts/install.py --prefix /user/local/bin
```

If anything goes wrong with the build process and you need to start again, be sure to add the command `--clean`.

You can also request all dependencies to be installed locally:

```shell
$ ./scripts/install.py --dependencies
```

## CMake

If Python3 isn't available, the binaries can be installed manually with [CMake](https://cmake.org):

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git
$ cd octopus/build
$ cmake .. && make install
```

To install to root (e.g. `/usr/local/bin`) use the `-D` option:

```shell
$ cmake -DINSTALL_ROOT=ON ..
```

CMake will try to find a suitable compiler on your system, if you'd like you use a specific compiler use the `-D` option, for example:

```shell
$ cmake -D CMAKE_C_COMPILER=clang-4.0 -D CMAKE_CXX_COMPILER=clang++-4.0 ..
```

You can check installation was successful by executing the command:

```shell
$ octopus -h
```

## Docker

## Conda

Octopus is available [pre-built for Linux](https://anaconda.org/bioconda/octopus) as part of [Bioconda](https://bioconda.github.io/). To [install in an isolated environment](https://bioconda.github.io/#using-bioconda):

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p venv
    venv/bin/conda install -c conda-forge -c bioconda octopus
    venv/bin/octopus -h

A package will also be available for OSX once conda-forge and bioconda move to newer versions of gcc and boost.