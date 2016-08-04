![Octopus Logo](logo.png)

[![Build Status](https://travis-ci.com/dancooke/octopus.svg?token=U9L3a7MWio2P3XpPT3JV&branch=master)](https://travis-ci.com/dancooke/octopus)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

---

**Warning: this project is incomplete and untested - do not use it for production work.**

---

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Octopus explicitly stores allele phasing which allows haplotypes to be dynamically excluded and extended whilst retaining the existing phasing information. Primarily this means Octopus can jointly consider allele sets far exceeding the cardinality of other approaches. But perhaps more importantly, this allows *marginalisation* over posterior distributions in haplotype space at specific loci. In practise this means Octopus can achieve far greater allelic genotyping accuracy than other methods, but can also infer conditional or unconditional phase probabilities directly from genotype probability distributions. This allows Octopus to report consistent allele event level variant calls *and* independent phase information.

##Requirements
* A C++14 compiler with SSE2 support
* Git 2.5 or greater
* Boost 1.61 or greater
* htslib 1.31.1 or greater
* CMake 3.3 or greater
* Optional:
    * Python3 or greater

##Installation

Octopus can be built and installed on a wide range of operating systems including most Unix based systems (Linux, OS X) and Windows.

####*Installing with Homebrew (for OS X)*
The recommended method of installation for Max OS X is with the package manager [Homebrew](http://brew.sh)

```shell
$ brew tap science
$ brew install octopus
```

This will download the relevant files (including any dependancies) and install to `usr/local/bin`. Octopus can then be easily updated or removed

```shell
$ brew upgrade octopus
$ brew uninstall octopus
```

####*Quick installation with Python3*
Manually installing octopus first requires obtaining the binaries. First `cd` to a writable directory and execute

```shell
$ git clone https://github.com/dancooke/octopus.git
```

then use the `Python3` install helper

```shell
$ ./octopus/make.py
```

####*Installing with CMake*
If `Python3` isn't available, the binaries can be installed manually with [CMake](https://cmake.org)

```shell
$ git clone https://github.com/dancooke/octopus.git
$ cd octopus/build
$ cmake .. && make install
```

##Running Tests

Test the installation was successful by executing the command 

```shell
$ octopus --help
```

####*Running the tests with Python3*
```shell
$ ./test.py
```

####*Running tests with CMake*
```shell
$ cd build
$ cmake -DBUILD_TESTING=ON .. && make test
```

##Examples

####*Calling germline variants in an individual*
```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam
```

####*Joint variant calling*
```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam NA12891.bam NA12892.bam
```

####*Targeted calling*
```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam --regions 1 2:30,000,000- 3:10,000,000-20,000,000
```

####*Calling somatic variants*
```shell
$ octopus --caller somatic --reference hs37d5.fa --reads normal.bam tumour.bam --normal-sample NORMAL
```

####*Calling denovo mutations*
```shell
$ octopus --caller denovo --reference hs37d5.fa --reads NA12878.bam NA12891.bam NA12892.bam --maternal-sample NA12892 --paternal-sample NA12891
```

##Documentation

Complete user and developer documentation is available in the doc directory.

##Support

Please report any bugs or feature requests to the [octopus issue tracker](https://github.com/dancooke/octopus/issues).

##Contributing

Contributions are very welcome, but please first review the [contribution guidelines](CONTRIBUTING.md).

##Authors

Daniel Cooke & Gerton Lunter

##Citing

##License

Please refer to [LICENSE](LICENSE).
