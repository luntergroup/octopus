---
id: germline
title: Germline WGS
---

This case study looks at whole-genome germline variant calling in a single individual (NA12878).

### Prerequisites

- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [BWA](https://github.com/lh3/bwa)
- [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools)
- [Octopus](https://github.com/luntergroup/octopus) (with [random forests](https://github.com/luntergroup/octopus/wiki/Variant-filtering:-Random-Forest) installed).

### Download data files

First download raw reads from the Illumina platinum genomes project for individual NA12878:

```shell
$ mkdir ~/fastq && cd ~/fastq
$ wget https://storage.googleapis.com/genomics-public-data/platinum-genomes/fastq/ERR194147_1.fastq.gz
$ wget https://storage.googleapis.com/genomics-public-data/platinum-genomes/fastq/ERR194147_2.fastq.gz
```

Next download a copy of the human reference sequence. In this example we use GRCh37 plus a decoy contig (recommended). If you prefer to use GRCh38, be sure to get a copy without alternative contigs or patches (but with a decoy contig), such as the one available [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz).

```shell
$ mkdir ~/reference && cd ~/reference
$ wget ftp://ftp-$ trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && gzip -d hs37d5.fa.gz
```

To evaluate our calls we need a truth set. We use the Genome In a Bottle (GIAB) version 3.3.2 high confidence calls for NA12878 (HG001):

```shell
$ mkdir ~/vcf/giab && cd ~/vcf/giab
$ wget ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37
```

### Map reads to reference genome

First we need to index the reference sequence:

```shell
$ cd ~/reference
$ samtools faidx hs37d5.fa
$ bwa index hs37d5.fa
```

Then map our reads to the reference:

```shell
$ mkdir ~/bam 
$ bwa mem -t 15 -R "@RG\tID:1\tSM:NA12878\tLB:platinum\tPU:illumina" \
      ~/reference/hs37d5.fa \
      ~/fastq/ERR194147_1.fastq.gz ~/fastq/ERR194147_2.fastq.gz \
      | samtools view -bh | samtools sort -@ 5 -o ~/bam/NA12878.platinum.b37.bam -
$ samtools index ~/bam/NA12878.platinum.b37.bam
```

This should take a few hours.

### Call variants

We do not recommend pre-processing the raw BWA alignments (e.g. duplicate marking, or base quality score recalibration) as we do not find this provides consistent  improvements in accuracy, and tends to slow down calling as pre-processed reads files can be larger than the originals. As this is human data, the default arguments for octopus should work well. We restrict calling to the autosomes plus X as these are the only contigs present in the validation sets. We also request a 'legacy' VCF file to use for benchmarking (see section on octopus's default VCF format).

```shell
$ octopus \
     -R ~/reference/hs37d5.fa \
     -I ~/bam/NA12878.platinum.b37.bam \
     -T 1 to MT \
     --sequence-error-model PCR-FREE.HISEQ-2000 \
     --forest forests/germline.v0.6.0-beta.forest \
     -o ~/vcf/NA12878.platinum.b37.octopus.vcf.gz \
     --threads 20 \
     --legacy
```

This should complete in a few hours.

### Evaluate variant calls

Finally, we will evaluate our calls with RTG Tools `vcfeval`. This command requires the reference sequence to be preprocessed:

```shell
$ ~/tools/rtgtools/rtg format -o ~/reference/hs37d5_sdf ~/reference/hs37d5.fa
```

Then run vcfeval:

```shell
$ rtg vcfeval \
     -t ~/reference/hs37d5_sdf \
     -b ~/vcf/giab/HG001_GRCh37_truth.vcf.gz \
     --evaluation-regions ~/vcf/giab/HG001_GRCh37_hiconf.bed \
     -c ~/vcf/NA12878.platinum.b37.octopus.legacy.vcf.gz \
     -o ~/eval/NA12878.platinum.b37.octopus.eval \
     -f FORMAT.RFQUAL \
     --ref-overlap
```

We see the following results:

```shell
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    4.770            3680280        3703892       2042      10581     0.9994       0.9971     0.9983
     None            3681055        3704782       3297       9806     0.9991       0.9973     0.9982
```