![Octopus Logo](logo.png)

[![Build Status](https://travis-ci.org/luntergroup/octopus.svg?branch=develop)](https://travis-ci.org/luntergroup/octopus)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![Gitter](https://badges.gitter.im/octopus-caller/Lobby.svg)](https://gitter.im/octopus-caller/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

---

**Warning: this project is incomplete - it may be unstable and contain bugs.**

---

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Octopus explicitly stores allele phasing infomation which allows haplotypes to be dynamically excluded and extended. Primarily this means octopus can jointly consider allele sets far exceeding the cardinality of other approaches, but perhaps more importantly, it allows *marginalisation* over posterior distributions in haplotype space at specific loci. In practise this means octopus can achieve far greater allelic genotyping accuracy than other methods, but can also infer conditional or unconditional phase probabilities directly from genotype probability distributions. This allows octopus to report consistent allele event level variant calls *and* independent phase information.

##Requirements
* A C++14 compiler with SSE2 support
* A C++14 standard library implementation
* Git 2.5 or greater
* Boost 1.58 or greater
* htslib 1.31.1 or greater
* CMake 3.5 or greater
* Optional:
    * Python3 or greater
    
Warning: GCC 6.1.1 and below have bugs which affect octopus, the code may compile, but do not trust the results. GCC 6.2 should be safe. Clang 3.8 has been tested. Visual Studio likely won't compile as it is not C++14 feature complete.

####*Obtaining requirements on OS X*

On OS X, Clang is recommended. All requirements can be installed using the package manager [Homebrew](http://brew.sh/index.html):

```shell
$ brew update
$ brew install git
$ brew install --with-clang llvm
$ brew install boost
$ brew install cmake
$ brew tap homebrew/science # required for htslib
$ brew install htslib
$ brew install python3
```
Note if you already have any of these packages installed via Homebrew on your system the command will fail, but you can update to the latest version using `brew upgrade`.

####*Obtaining requirements on Ubuntu Xenial*

On Ubuntu, Clang 3.8 is recommended as GCC 6.2 has a [bug](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77550) which slows down octopus. To install the requirements use:

```shell
$ sudo apt-get update && sudo apt-get upgrade
$ sudo apt-get install clang-3.8
$ sudo apt-get install libstdc++6
$ sudo apt-get install libboost-all-dev
$ sudo apt-get install cmake
$ sudo apt-get install git-all
$ sudo apt-get install python3
```

Only htslib 1.2.1 is currently available via `apt-get` so this will need to be installed manually:

```shell
$ sudo apt-get install autoconf
$ git clone https://github.com/samtools/htslib.git
$ autoheader
$ autoconf
$ ./configure
$ make && sudo make install
```

##Installation

Octopus can be built and installed on a wide range of operating systems including most Unix based systems (Linux, OS X) and Windows (once MSVC is C++14 feature complete).

####*Quick installation with Python3*

Manually installing octopus requires obtaining a copy the binaries. In the command line, direct to an appropriate install directory and execute:

```shell
$ git clone https://github.com/luntergroup/octopus.git
```

then *if your default compiler is Clang 3.8 or GCC 6.2* use the Python3 install helper:

```shell
$ ./octopus/install.py
```

otherwise explicitly specify the compiler to use:

```shell
$ ./octopus/install.py --compiler /path/to/cpp/compiler # or just the compiler name if on your PATH
```

For example, if the requirement instructions above were used:

```shell
$ ./octopus/install.py --compiler clang++-3.8
```

By default this installs to `/bin` relative to where you installed octopus. To install to a root directory (e.g. `/usr/local/bin`) use:

```shell
$ ./octopus/install.py --root
```

this may prompt you to enter a `sudo` password.

If anything goes wrong with the build process and you need to start again, be sure to add the command `--clean`!

####*Installing with CMake*

If Python3 isn't available, the binaries can be installed manually with [CMake](https://cmake.org):

```shell
$ git clone https://github.com/luntergroup/octopus.git
$ cd octopus/build
$ cmake ..
$ make install
```

By default this installs to the `/bin` directory where octopus was installed. To install to root (e.g. `/usr/local/bin`) use the `-D` option:

```shell
$ cmake -DINSTALL_ROOT=ON ..
```

CMake will try to find a suitable compiler on your system, if you'd like you use a specific compiler use the `-D` option, for example:

```shell
$ cmake -D CMAKE_C_COMPILER=clang-3.8 -D CMAKE_CXX_COMPILER=clang++-3.8 ..
```

You can check installation was successful by executing the command:

```shell
$ octopus --help
```

##Running Tests

Octopus comes packaged with unit, regression, and benchmark testing. The unit tests are self-contained whilst the regression and benchmark tests require external data sources.

####*Running the tests with Python3*

```shell
$ test/install.py
```

####*Running tests with CMake*

```shell
$ cd build
$ cmake -DBUILD_TESTING=ON .. && make test
```

##Examples

Here are some common use-cases to get started. These examples are by no means exhaustive, please consult the documentation for explanations of all options, algorithms, and further examples.

Note by default octopus will output all calls in VCF format to standard output, in order to write calls to a file (`.vcf`, `.vcf.gz`, and `.bcf` are supported), use the command line option `--output` (`-o`).

####*Calling germline variants in an individual*

This is the simplest case, if the file `NA12878.bam` contains a single sample, octopus will default to its individual calling model:

```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam
```

or less verbosely:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam
```

Note octopus automatically detects the samples contained in the input read file and will jointly call (see below) all samples present by default, to restrict calling to a single sample in this case it is required to explicitly specify which sample to use:

```shell
$ octopus -R hs37d5.fa -I multi-sample.bam -S NA12878
```

####*Targeted calling*

By default, octopus will call all possible regions (as specified in the reference FASTA). In order to select a set of target regions, use the `--regions` (`-T`) option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -T 1 2:30,000,000- 3:10,000,000-20,000,000
```

Or conversely a set of regions to *exclude* can be given with `--skip-regions` (`-K`):

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -K 1 2:30,000,000- 3:10,000,000-20,000,000
```

####*Calling de novo mutations in a trio*

To call germline and de novo mutations in a trio, just specify maternal (`--maternal-sample`; `-M`) and paternal (`--paternal-sample`; `-F`) samples:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -M NA12892 -F NA12891
```

####*Calling somatic mutations in tumours*

By default octopus will use either the individual or population models. To invoke another calling model, either the `--caller` (`-C`) option must be set, or another option which allows octopus to deduce the required caller must be supplied.

In the case of somatic variant, the somatic caller is automatically invoked if the option `--normal-sample` (`-N`) is set:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam --normal-sample NORMAL
```

This will call somatic and germline variants in a single individual.
 
It is also possible to call multiple tumours from the same individual jointly:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumourA.bam tumourB --normal-sample NORMAL
```

Octopus will then emit separate genotype calls for each sample.

If a normal sample is not present the cancer model must be invoked explicitly:

```shell
$ octopus -R hs37d5.fa -I tumour1.bam tumour2.bam -C cancer
```

Note without a normal sample classification power is significantly reduced.

####*Joint variant calling*

Octopus uses different calling models for populations and individuals. Briefly, the individual model is exact whilst the population model uses approximations. However, it is recommended to use the population model to call *germline variants* in multiple samples from the same population as the model can leverage information between individuals:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam
```

####*Multithreaded calling*

All the octopus callers can be run using multiple cores, activated with the `--threads` command line option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads
```

This will let octopus automatically decide how many threads to use, and is the recommended approach as octopus can dynamically juggle thread usage at an algorithm level. However, a strict upper limit on the number of threads can also be used:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads=4
```

####*Fast calling*

By default, octopus is geared towards more accurate variant calling which requires the use of complex (slow) algorithms. If speed is a concern for you, then many of these features can be disabled to get very fast runtimes:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --fast
```

Note this does not turn on multithreading or increase buffer sizes.

##Documentation

Complete user and developer documentation is available in the doc directory.

##Support

Please report any bugs or feature requests to the [octopus issue tracker](https://github.com/dancooke/octopus/issues).

##Contributing

Contributions are very welcome, but please first review the [contribution guidelines](CONTRIBUTING.md).

##Authors

Daniel Cooke and Gerton Lunter

##Citing

TBA

##License

Please refer to [LICENSE](LICENSE).
