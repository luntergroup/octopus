![Octopus Logo](logo.png)

[![Website](https://img.shields.io/website?url=https%3A%2F%2Fluntergroup.github.io%2Foctopus%2F)](https://luntergroup.github.io/octopus/)
[![Build Status](https://travis-ci.org/luntergroup/octopus.svg?branch=master)](https://travis-ci.org/luntergroup/octopus)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![Gitter](https://badges.gitter.im/octopus-caller/Lobby.svg)](https://gitter.im/octopus-caller/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/luntergroup/octopus)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/octopus/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/dancooke/octopus?label=docker)](https://hub.docker.com/r/dancooke/octopus)

Octopus is a mapping-based variant caller that implements several calling models within a unified haplotype-aware framework. Octopus takes inspiration from particle filtering by constructing a tree of haplotypes and dynamically pruning and extending the tree based on haplotype posterior probabilities in a sequential manner. This allows octopus to implicitly consider all possible haplotypes at a given loci in reasonable time.

There are currently six calling models available:

- [individual](https://luntergroup.github.io/octopus/docs/guides/models/individual): call germline variants in a single healthy individual.
- [population](https://luntergroup.github.io/octopus/docs/guides/models/population): jointly call germline variants in small cohorts.
- [trio](https://luntergroup.github.io/octopus/docs/guides/models/trio): call germline and _de novo_ mutations in a parent-offspring trio.
- [cancer](https://luntergroup.github.io/octopus/docs/guides/models/cancer): call germline and somatic mutations tumour samples.
- [polyclone](https://luntergroup.github.io/octopus/docs/guides/models/polyclone): call variants in samples with an unknown mixture of haploid clones, such a bacteria or viral samples.
- [cell](https://luntergroup.github.io/octopus/docs/guides/models/cell): call variants in a set of single cell samples from the same individual.

Octopus calls SNVs, small-medium sized indels, and small complex rearrangements in [VCF 4.3](https://luntergroup.github.io/octopus/docs/guides/advanced/vcf).

## Quick start

Install Octopus:

```shell
$ git clone https://github.com/luntergroup/octopus.git
$ octopus/scripts/install.py --dependencies --forests
$ echo 'export PATH='$(pwd)'/octopus/bin:$PATH' >> ~/.bash_profile
$ source ~/.bash_profile
```

Call some variants:

```shell
$ FOREST="$(pwd)/octopus/resources/forests/germline.v0.7.4.forest"
$ octopus -R hs37d5.fa -I NA12878.bam -T 1 to MT -o NA12878.octopus.vcf.gz --forest $FOREST --threads 8
```

## Documentation

Documentation is hosted on [GitHub pages](https://luntergroup.github.io/octopus/).

## Support

Please report any bugs or feature requests to the [octopus issue tracker](https://github.com/luntergroup/octopus/issues). General chat is hosted on [Gitter](https://gitter.im/octopus-caller/Lobby).

## Contributing

Contributions are very welcome! Please consult the [contribution guidelines](CONTRIBUTING.md).

## Authors

Daniel Cooke and Gerton Lunter

## Citing

See [publications](https://luntergroup.github.io/octopus/docs/publications) associated with Octopus.

## License

Octopus is distributed under the [MIT LICENSE](LICENSE).
