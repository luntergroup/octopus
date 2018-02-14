![Octopus Logo](logo.png)

[![Build Status](https://travis-ci.org/luntergroup/octopus.svg?branch=master)](https://travis-ci.org/luntergroup/octopus)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![Gitter](https://badges.gitter.im/octopus-caller/Lobby.svg)](https://gitter.im/octopus-caller/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/octopus/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

---

**Warning: this project is incomplete - it may be unstable and contain bugs.**

---

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Octopus takes inspiration from particle filtering by constructing a tree of haplotypes and dynamically pruning and extending the tree based on haplotype posterior probabilities in a sequential manner. This allows octopus to implicitly candider all possible haplotypes at a given loci in reasonable time.

There are currently three calling models implemented:

- An individual model for calling **germline variants** in a single healthy individual.
- A tumour model for calling germline variants and **somatic mutations** in one or more tumour samples from a single individual.
- A trio model for calling germline variants and *de novo* mutations in a parent-offspring trio.

Octopus is currently able to call SNVs, small-medium sized indels, small complex rearrangements and micro-inversions.

We hope to implement more calling models in the future, including, but not limited to:

- A population model for calling germline variants from multiple healthy individuals within a population.
- A pedigree model for calling germline and *de novo* mutations in multiple healthy individuals from a known pedigree.
- A haploid clonal model for calling polymorphisms in a sample of mixed bacteria isolates.

## Requirements
* A C++14 compiler with SSE2 support
* A C++14 standard library implementation
* Git 2.5 or greater
* Boost 1.65 or greater
* htslib 1.4 or greater
* CMake 3.9 or greater
* Optional:
    * Python3 or greater

#### *Obtaining requirements on OS X*

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

#### *Obtaining requirements on Ubuntu*

Depending on your Ubuntu distribution, some requirements can be installed with `apt-get`. It may be preferable to use GCC as this will simplify installing Boost:

```shell
$ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
$ sudo apt-get update && sudo apt-get upgrade
$ sudo apt-get install gcc-7
$ sudo apt-get install git-all
$ sudo apt-get install python3
```

The other packages will need to be installed manually:

- CMake installation instructions are given [here](https://askubuntu.com/a/865294).
- Htslib installation instructions are given [here](https://github.com/samtools/htslib). Note you may need to install `autoconf` (`sudo apt-get install autoconf`).
- Instructions on installing Boost are given [here](https://stackoverflow.com/a/24086375/2970186).

These instructions are replicated in the [user documentation](https://github.com/luntergroup/octopus/blob/develop/doc/manuals/user/octopus-user-manual.pdf) (Appendix).

## Installation

Octopus can be built and installed on a wide range of operating systems including most Unix based systems (Linux, OS X) and Windows (once MSVC is C++14 feature complete).

#### Conda package

Octopus is available [pre-built for Linux](https://anaconda.org/bioconda/octopus) as part of [Bioconda](https://bioconda.github.io/). To [install in an isolated environment](https://bioconda.github.io/#using-bioconda):

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p venv
    venv/bin/conda install -c conda-forge -c bioconda octopus
    venv/bin/octopus -h

A package will also be available for OSX once conda-forge and bioconda move to newer versions of gcc and boost.

#### *Quick installation with Python3*

First clone the git repository in your preferred directory:

```shell
$ git clone -b master https://github.com/luntergroup/octopus.git && cd octopus
```

The easiest way to install octopus from source is with the Python3 install script. If your default compiler satisfies the minimum requirements just execute:

```shell
$ ./install.py
```

otherwise explicitly specify the compiler to use:

```shell
$ ./install.py --cxx_compiler /path/to/cpp/compiler # or just the compiler name if on your PATH
```

For example, if the requirement instructions above were used:

```shell
$ ./install.py --cxx_compiler clang++-4.0
```

On some systems, you may also need to specify a C compiler which is the same version as your C++ compiler, otherwise you'll get lots of link errors. This can be done with the `--c_compiler` option:

```shell
$ ./install.py -cxx g++-7 -c gcc-7 
```

By default this installs to `/bin` relative to where you installed octopus. To install to a root directory (e.g. `/usr/local/bin`) use:

```shell
$ ./install.py --root
```

If anything goes wrong with the build process and you need to start again, be sure to add the command `--clean`.

#### *Installing with CMake*

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

## Running Tests

Octopus currently has limited unit tests (more are planned!). To install and run them, use the Python3 install script in the `test` directory:

```shell
$ test/install.py
```

## Examples

Here are some common use-cases to get started. These examples are by no means exhaustive, please consult the documentation for explanations of all options, algorithms, and further examples. For a more in depth example, refer to the [whole genome germline calling case study](https://github.com/luntergroup/octopus/blob/master/doc/octopus_wgs_case_study.md).

Note by default octopus will output all calls in VCF format to standard output, in order to write calls to a file (`.vcf`, `.vcf.gz`, and `.bcf` are supported), use the command line option `--output` (`-o`).

#### *Calling germline variants in an individual*

This is the simplest case, if the file `NA12878.bam` contains a single sample, octopus will default to its individual calling model:

```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam
```

or less verbosely:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam
```

By default, octopus automatically detects and calls all samples contained in the input read files. To call a subset of these samples, use the `--samples` (`-S`) option:

```shell
$ octopus -R hs37d5.fa -I multi-sample.bam -S NA12878
```

#### *Targeted calling*

By default, octopus will call all regions specified in the reference index. In order to restrict calling to a subset of regions, either provide a list of zero-indexed regions in the format `chr:start-end` (`--regions`; `-T`), or a file containing a list of regions in either standard format or BED format (`--regions-file`; `-t`):

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -T 1 2:30,000,000- 3:10,000,000-20,000,000
$ octopus -R hs37d5.fa -I NA12878.bam -t regions.bed
```

Conversely a set of regions to *exclude* can be given explictely (`--skip-regions`;`-K`), or with a file (`--skip-regions-file`; `-k`):

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -K 1 2:30,000,000- 3:10,000,000-20,000,000
$ octopus -R hs37d5.fa -I NA12878.bam -k skip-regions.bed
```

#### *Calling de novo mutations in a trio*

To call germline and de novo mutations in a trio, either specify both maternal (`--maternal-sample`; `-M`) and paternal (`--paternal-sample`; `-F`) samples:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -M NA12892 -F NA12891
```

or provide a PED file which defines the trio:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam --pedigree ceu_trio.ped
```

#### *Calling somatic mutations in tumours*

To call germline and somatic mutations in a paired tumour-normal sample, just specify which sample is the normal (`--normal-sample`; `-N`):

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam --normal-sample NORMAL
```
 
It is also possible to genotype multiple tumours from the same individual jointly:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumourA.bam tumourB.bam --normal-sample NORMAL
```

If a normal sample is not present the cancer calling model must be invoked explicitly:

```shell
$ octopus -R hs37d5.fa -I tumour1.bam tumour2.bam -C cancer
```

Note however, that without a normal sample, somatic mutation classification power is significantly reduced.

#### *Joint variant calling (experimental)*

Multiple samples from the same population, without pedigree information, can be called jointly:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam
```

Joint calling samples may increase calling power, especially for low coverage sequencing.

#### *HLA genotyping*

To call phased HLA genotypes, increase the default phase level:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -t hla-regions.bed -l aggressive
```

#### *Multithreaded calling*

Octopus has built in multithreading capabilities, just add the `--threads` command:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads
```

This will let octopus automatically decide how many threads to use, and is the recommended approach as octopus can dynamically juggle thread usage at an algorithm level. However, a strict upper limit on the number of threads can also be used:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads 4
```

#### *Fast calling*

By default, octopus is geared towards more accurate variant calling which requires the use of complex (slow) algorithms. However, to achieve faster runtimes (at the cost of decreased calling accuracy) many of these features can be disabled. There are two helper commands that setup octopus for faster variant calling, `--fast` and `--very-fast`, e.g.:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --fast
```

Note this does not turn on multithreading or increase buffer sizes.

## Output format

Octopus outputs variants using a simple but rich VCF format (see [user documentation](https://github.com/luntergroup/octopus/blob/develop/doc/manuals/user/octopus-user-manual.pdf) for full details). For example, two overlapping deletions are represented like:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTAT	A,*	.	.	.	GT	1|2
1	102738191	.	ATTATTTATTTAT	A	.	.	.	GT	.|1
```

in contrast to how such a site would usually be represented, either:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTAT	A	.	.	.	GT	1/0
1	102738191	.	ATTATTTATTTAT	A	.	.	.	GT	1/0
```

which is inconsistent as the reference is deduced in each record, or:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTATTTAT	ATTAT,A	.	.	.	GT	1/2
```

which is at least consistent, but rapidly becomes unmanageable as the length and number of overlapping variants increases.

Octopus's representation is both succinct and consistent. The `*` allele denotes an upstream deletion, while the `.` in the genotype of the second record indicates the allele is missing due to a previous event. As the records are phased, the called haplotypes can be unambiguously reconstructed when the VCF file is read sequentially.

However, some existing tools will not recognise this format. For example, RTG Tools does not fully support this representation. Therefore, octopus has an option to also produce calls using a more typical VCF format (like the first of the two examples). To request this, use the `--legacy` command line option. This option is only available when outputting calls to a file (i.e. not `stdout`).   

## Documentation

Complete [user](https://github.com/luntergroup/octopus/blob/develop/doc/manuals/user/octopus-user-manual.pdf) and [developer](https://github.com/luntergroup/octopus/blob/develop/doc/manuals/dev/octopus-dev-manual.pdf) documentation is available in the doc directory.

## Support

Please report any bugs or feature requests to the [octopus issue tracker](https://github.com/luntergroup/octopus/issues).

## Contributing

Contributions are very welcome, but please first review the [contribution guidelines](CONTRIBUTING.md).

## Authors

Daniel Cooke and Gerton Lunter

## Citing

TBA

## License

Please refer to [LICENSE](LICENSE).
