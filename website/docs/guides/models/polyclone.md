---
id: polyclone
title: Polyclone
---

The `polyclone` calling model is designed for calling variation in a pooled sample of haploid genomes where the number and mixture composition of genomes is unknown. The model will attempt to deconvolute haplotypes and call the number of local haplotypes, in addition to variants on called haplotypes. Use cases include bacterial or viral metagenomics, mixed infection analysis, and mitochondria sequencing.

## Usage

The polyclone calling model is invoked by speciying the `polyclone` [`--caller`](cli.md#--caller).

```shell
$ octopus -R reference.fa -I reads.bam -C polyclone
```

## VCF output



## Performance considerations

The most important parameter in this model is `--max-clones` which determines the maximum number of haplotypes that octopus will try to use to fit the data. This is a bit like setting the 'ploidy' for the sample. Higher values of `--max-clones` may lead to longer runtimes and more memory usage - we do not recommend values above 10.