---
id: individual
title: individual
---

The `individual` calling model is used to call germline variants in a single sample of known ploidy. It is the simplest model.

## Usage

If the file `NA12878.bam` contains a single sample, to call germline variants in all regions use:

```shell
$ octopus --reference hs37d5.fa --reads NA12878.bam
```

or less verbosely:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam
```

By default, octopus automatically detects and calls all samples contained in the input read files. If your BAM file contains multiple samples, and you want to call just one of these, use the `--samples` (`-S`) option:

```shell
$ octopus -R hs37d5.fa -I multi-sample.bam -S NA12878
```

## Setting ploidy

Octopus is setup to call diploid samples by default. If your sample is not diploid you can set the ploidy with the `--organism-ploidy` (`-P`) option:

```shell
$ octopus -R hs37d5.fa -I haploid.bam -P 1
```

You can also set contig specific policies with the `--contig-ploidies` (`-p`) option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -p Y=1
```

Octopus automatically sets contigs `Y` and `chrY` to ploidy 1.

_Performance note_: Larger ploidies require greater computational resources. In general, ploidies above 4 are currently intractable. If you wish to call tetraploid sample, you may find that you need to tweak other performance related parameters to get reasonable run times.