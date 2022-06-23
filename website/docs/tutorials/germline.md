---
id: germline
title: Germline WGS
---

This case study considers whole-genome germline variant calling in an individual. We will use the syntheic-diploid sample [CHM1-CHM13](https://www.nature.com/articles/s41592-018-0054-7).

This tutorial is [included](../../static/snakemake/germline.smk) as a [Snakemake](https://snakemake.readthedocs.io/en/stable/#) workflow. Run it with

```shell
$ snakemake --snakefile germline.smk -j 16 --use-singularity
```

## Prerequisites

- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [BWA](https://github.com/lh3/bwa)
- [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools)
- [Octopus](https://github.com/luntergroup/octopus) (with [random forests](guides/../../guides/filtering/forest.md) installed).

This tutorial assumes a directory structure like

```
.
└───data
│   └───references
│   └───reads
│   │   └───raw
│   │   └───mapped
│   └───truth
└───results
    └───calls
    └───eval
```

We'll go ahead and create this upfront:

```shell
$ mkdir -p data/{references,truth} data/reads/{raw,mapped}
$ mkdir -p results/{calls,eval}
```

## Download data

First, download WGS CHM1-CHM13 PCR-free Illumina HiSeq-X10 reads ([PRJEB13208](https://www.ebi.ac.uk/ena/browser/view/PRJEB13208)) from EBI:

```shell
$ curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_1.fastq.gz | gzip > data/reads/raw/CHM1-CHM13_1.fastq.gz
$ curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_2.fastq.gz | gzip > data/reads/raw/CHM1-CHM13_2.fastq.gz
```

Next, we need a reference genome. In this example we'll use GRCh38 with ALT and decoys contigs (hs38DH). The simplest way to get this is with `bwakit`:

```shell
$ run-gen-ref hs38DH
$ mv hs38DH.fa* data/references
```

We'll also need the truth variants to evaluate our calls, which we can get from the [CHM-eval GitHub](https://github.com/lh3/CHM-eval):

```shell
$ curl -L https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar | tar xf - > data/truth/syndip
```

## Map reads

First, we need to index the reference genome for alignment:

```shell
$ samtools faidx data/references/hs38DH.fa
$ bwa index data/references/hs38DH.fa
```

Next, we map the raw reads to the reference genome with `bwa mem`:

```shell
$ bwa mem -t 16 \
     -R "@RG\tID:0\tSM:CHM1-CHM1378\tLB:Broad\tPU:Illumina" \
     data/reference/hs38DH.fa \
     data/reads/raw/CHM1-CHM13_1.fastq.gz data/reads/raw/CHM1-CHM13_2.fastq.gz | \
     samtools view -bh | \
     samtools sort -@ 4 -o data/reads/mapped/CHM1-CHM13.hs38DH.bwa-mem.bam -
$ samtools index data/reads/mapped/CHM1-CHM13.hs38DH.bwa-mem.bam
```

This should take a few hours.

## Call variants

Now we can call variants with `octopus`. Since this is single-sample germline calling, we'll be using the [individual](../guides/models/individual.md) calling model. We'll also set the [sequence error model](../guides/errorModels.md) to reflect the PCR-free library design of this sample and Illlumina X10 sequencer, and use [random forest filtering](../guides/filtering/forest.md). Finally, we [restrict calling](../guides/advanced/targeted.md) to the primary chromosomes and use the built-in [multithreading](../guides/advanced/threading.md):

```shell
$ octopus \
     -R data/reference/hs38DH.fa \
     -I data/reads/mapped/CHM1-CHM13.hs38DH.bwa-mem.bam \
     -T chr1 to chrM \
     --sequence-error-model PCRF.X10 \
     --forest resources/forests/germline.v0.8.0.forest \
     -o results/calls/CHM1-CHM13.hs38DH.bwa-mem.octopus.vcf.gz \
     --threads 16
```

This should complete in a few hours.

## Evaluate variants

Finally, we will evaluate our calls with RTG Tools `vcfeval`. This tool requires the reference sequence to be preprocessed:

```shell
$ rtg format -o ~/reference/hs37d5.sdf ~/reference/hs37d5.fa
```

Then evaluate the calls with `vcfeval`:

```shell
$ rtg vcfeval \
     -t data/reference/hs38DH.sdf \
     -b data/truth/syndip/full.38.vcf.gz \
     --evaluation-regions data/truth/syndip/full.38.bed.gz \
     -c results/calls/CHM1-CHM13.hs38DH.bwa-mem.octopus.vcf.gz \
     -o results/eval/CHM1-CHM13.hs38DH.bwa-mem.octopus.pass.vcfeval \
     -f RFGQ \
     -m annotate \
     --ref-overlap
```

We should see the following results:

```shell
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
```