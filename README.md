![Octopus Logo](logo.png)

[![Build Status](https://travis-ci.com/dancooke/octopus.svg?token=U9L3a7MWio2P3XpPT3JV&branch=master)](https://travis-ci.com/dancooke/octopus)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

---

**Warning: this project is incomplete and untested - do not use it for production work.**

---

**Octopus** is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Unlike other variant callers, Octopus explicitly stores allele phasing which allows haplotypes to be dynamically excluded and extended whilst retaining the existing phasing information. Primarily this means Octopus can jointly consider allele sets far exceeding the cardinality of other approaches. But perhaps more importantly, this allows *marginalisation* over posterior distributions in haplotype space at specific loci. In practise this means Octopus can achieve far greater allelic genotyping accuracy than other methods, but can also infer conditional or unconditional phase probabilities directly from genotype probability distributions. This allows Octopus to report consistent allele event level variant calls *and* independent phase information.

##Requirements
* A C++14 compiler
* Git 2.5 or greater
* Boost 1.61 or greater
* htslib 1.31.1 or greater
* CMake 3.3 or greater
* Optional:
    * Python3 or greater

##Installation

####*Installing with Homebrew (for MacOSX)*
```shell
$ brew install octopus
```

####*Quick installation with Python3*
```shell
$ ./make.py
```

####*Installing with CMake*
```shell
$ cd build
$ cmake .. && make install
```

Test the installation was successful by executing the command 

```shell
$ octopus --help
```

##Running tests

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

Daniel Cooke

##Citing



##License

Please refer to [LICENSE](LICENSE).
