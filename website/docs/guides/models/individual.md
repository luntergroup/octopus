---
id: individual
title: Individual
---

The `individual` calling model is used to call germline variants in a single sample with known ploidy. It is the simplest model Octopus offers.

## Basic usage

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

## Setting the ploidy

Octopus assume diploid samples by default. If your sample is not diploid you can set the ploidy with the `--organism-ploidy` (`-P`) option:

```shell
$ octopus -R hs37d5.fa -I haploid.bam -P 1
```

You can also set contig specific policies with the `--contig-ploidies` (`-p`) option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -p Y=1
```

:::important

There are `binom(|haplotypes| + ploidy - 1, ploidy)` genotypes for a given set of `haplotypes`, a number that rapidly becomes intractable for `ploidy > 2`. For polyploid samples, the number of genotypes considered is therefore limited to `--max-genotypes` using heuristics. The calling and runtime performance of the method, particularly phasing accuracy, can therefore be strongly dependent on `--max-genotypes` for polyploids.

:::