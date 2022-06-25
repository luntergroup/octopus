---
id: mitochondria
title: Mitochondria
---

This case study considers Mitochondria variant calling. We will call and benchmark variants in the mixture sample described in [Fazzini et al.](https://www.mdpi.com/1422-0067/22/2/935).

This tutorial is [included](../../static/snakemake/mitochondria.smk) as a [Snakemake](https://snakemake.readthedocs.io/en/stable/#) workflow. Run it with

```shell
$ snakemake --snakefile mitochondria.smk -j 16 --use-singularity
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

First, download PCR-amplified mitochondria Illumina MiSeq reads:

```shell
$ curl -o data/reads/raw/H1_U5.M4-Herk_S41.R1.fastq.gz https://zenodo.org/record/3991749/files/M4-Herk_S41_R1_001.fastq.gz
$ curl -o data/reads/raw/H1_U5.M4-Herk_S41.R2.fastq.gz https://zenodo.org/record/3991749/files/M4-Herk_S41_R2_001.fastq.gz
```

Next, we need a reference genome. For this tutorial we'll use GRCh38 with ALT and decoys contigs (hs38DH), which includes the revised Cambridge Reference Sequence (rCRS). The simplest way to get this is with `bwakit`:

```shell
$ run-gen-ref hs38DH
$ mv hs38DH.fa* data/references
```

We'll also need the truth variants to evaluate our calls, Fazzini et al. provide these in [Supplementary Table S7](https://www.mdpi.com/1422-0067/22/2/935/s1?version=1611043746). However, since we'll need a VCF for benchmarking, we have manually translated these variants into VCF, which can be downloaded using:

```shell
$ curl -o data/truth/H1_U5.vcf.gz https://github.com/luntergroup/octopus/blob/develop/website/static/assets/H1_U5.vcf.gz
$ curl -o data/truth/H1_U5.vcf.gz.tbi https://github.com/luntergroup/octopus/blob/develop/website/static/assets/H1_U5.vcf.gz.tbi
$ curl -o data/truth/H1_U5.highconf.bed https://github.com/luntergroup/octopus/blob/develop/website/static/assets/H1_U5.highconf.bed
```

## Map reads

First, we need to index the reference genome for alignment:

```shell
$ samtools faidx data/references/hs38DH.fa
$ bwa index data/references/hs38DH.fa
```

Next, we map the raw reads to the reference genome with `bwa mem`:

```shell
$ bwa mem \
     -R "@RG\tID:S41\tSM:H1_U5\tLB:M4\tPU:Illumina" \
     data/reference/hs38DH.fa \
     data/reads/raw/H1_U5.M4-Herk_S41.R1.fastq.gz data/reads/raw/H1_U5.M4-Herk_S41.R2.fastq.gz | \
     samtools view -bh | \
     samtools sort -o data/reads/mapped/H1_U5.M4-Herk_S41.hs38DH.bam -
$ samtools index data/reads/mapped/H1_U5.M4-Herk_S41.hs38DH.bam
```

This should complete in less than a minute.

## Call variants

Now we can call variants with `octopus`. Since we're aiming to call homoplasmic and heteroplasmic mitochondria variants, we'll be using the [polyclone](../guides/models/polyclone.md) calling model. We'll use the provided mitochondria [config](../guides/advanced/configs.md) which sets several options. We'll also set the [sequence error model](../guides/errorModels.md) to reflect the PCR amplified library design of this sample. Finally, we [restrict calling](../guides/advanced/targeted.md) to the mitochondria reference genome and use built-in [multithreading](../guides/advanced/threading.md):

```shell
$ octopus \
     -R data/reference/hs38DH.fa \
     -I data/reads/mapped/H1_U5.M4-Herk_S41.hs38DH.bam \
     -T chrM \
     --config /opt/octopus/resources/configs/mitochondria.config \
     --sequence-error-model PCR \
     -o results/calls/H1_U5.M4-PCR-Herk_S14.hs38DH.octopus.vcf.gz \
     --threads 8
```

This should complete in ~30 minutes.

## Evaluate variants

Finally, we will evaluate our calls with RTG Tools `vcfeval`. This tool requires the reference sequence to be preprocessed:

```shell
$ rtg format -o ~/reference/hs38DH.sdf ~/reference/hs38DH.fa
```

Then evaluate the calls with `vcfeval`:

```shell
$ rtg vcfeval \
     -t data/reference/hs38DH.sdf \
     -b data/truth/H1_U5.vcf.gz \
     --evaluation-regions data/truth/H1_U5.highconf.bed \
     -c results/calls/H1_U5.M4-PCR-Herk_S14.hs38DH.octopus.vcf.gz \
     -o results/eval/H1_U5.M4-PCR-Herk_S14.hs38DH.octopus.pass.alt.vcfeval \
     --sample ALT \
     -m annotate \
     --ref-overlap
```

We should see the following results:

```shell
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
     None                 38             38          1          4     0.9744       0.9048     0.9383
```
