---
id: trio
title: trio
---

The `trio` calling model is for calling germline variation and _de novo_ mutations in parent-offspring trios.

## Usage

To call germline and de novo mutations in a trio, either specify both maternal (`--maternal-sample`; `-M`) and paternal (`--paternal-sample`; `-F`) samples:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -M NA12892 -F NA12891
```

or provide a PED file which defines the trio:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam --pedigree ceu_trio.ped
```

## Setting ploidies

The trio calling model is currently only fully supported for diploid samples. You can set sample contig ploidies with the `--contig-ploidies` (`-p`) option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -M NA12892 -F NA12891 -p NA12891:X=1
``` 

## *De novo* mutations only

To call only `DENOVO` mutations, just add the `--denovos-only` command:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam --ped ceu_trio.ped --denovos-only
```

Note this only works for filtered calls, and should only be used if you do not plan on re-filtering the calls, since this will not be possible once the germline mutations are removed.

## Performance consideration

Like for the population model, `--max-joint-genotypes` has a large baring on overall accuracy and runtime. Increasing this parameter can improve sensitivity, but will also result in longer runtimes. The default value is set with accuracy in mind.