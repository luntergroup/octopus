# Whole genome germline case study

Here we will work through a real whole genome calling case study, from FASTQ to VCF evaluation. In addition to octopus, we make use of the following software tools:

* [samtools](http://samtools.sourceforge.net) (version 1.6)
* [BWA](http://bio-bwa.sourceforge.net) (versio 0.7.17-r1188)
* [RTG Tools](https://www.realtimegenomics.com/products/rtg-tools) (version 3.8.4)

## Download data files

First download raw reads from the Illumina platinum genomes project for individual NA12878:

```
mkdir ~/data/fastq && cd ~/data/fastq
wget https://storage.googleapis.com/genomics-public-data/platinum-genomes/fastq/ERR194147_1.fastq.gz
wget https://storage.googleapis.com/genomics-public-data/platinum-genomes/fastq/ERR194147_2.fastq.gz
```

Next download a copy of the human reference sequence. In this example we use GRCh37 plus a decoy contig (recommended). If you prefer to use GRCh38, be sure to get a copy without alternative contigs or patches (but with a decoy contig), such as the one available [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz).

```
mkdir ~/data/reference && cd ~/data/reference
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz && gzip -d hs37d5.fa.gz
```

To evaluate our calls we need a truth set. We use the Genome In a Bottle (GIAB) version 3.3.2 high confidence calls for NA12878 (HG001):

```
mkdir ~/data/vcf/giab && cd ~/data/vcf/giab
wget ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37
```

## Map reads to reference genome

First we need to index the reference sequence:

```
cd ~/data/reference
samtools faidx hs37d5.fa
bwa index hs37d5.fa
```

Then map our reads to the reference:

```
mkdir ~/data/bam 
bwa mem -t 15 -R "@RG\tID:NA12878\tSM:NA12878\tLB:platinum\tPU:illumina" \
    ~/data/reference/hs37d5.fa \
    ~/data/fastq/ERR194147_1.fastq.gz ~/data/fastq/ERR194147_2.fastq.gz \
    | samtools view -bh > ~/data/bam/NA12878.platinum.b37.unsorted.bam
samtools sort -@ 15 -o ~/data/bam/NA12878.platinum.b37.bam ~/data/bam/NA12878.platinum.b37.unsorted.bam
samtools index ~/data/bam/NA12878.platinum.b37.bam
rm ~/data/bam/NA12878.platinum.b37.unsorted.bam
```

## Call variants

We do not recommend pre-processing the raw BWA alignments (e.g. duplicate marking, or base quality score recalibration) as we do not find this provides consistent  improvements in accuracy, and tends to slow down calling as pre-processed reads files are often considerably larger than the originals. As this is human data, the default arguments for octopus should work well. We restrict calling to the autosomes plus X as these are the only contigs present in the validation sets. We also request a 'legacy' VCF file to use for benchmarking (see section on octopus's default VCF format).

```
octopus -R ~/data/reference/hs37d5.fa \
    -I ~/data/bam/NA12878.platinum.b37.bam \
    -o ~/data/vcf/NA12878.platinum.b37.octopus.vcf.gz \
    -T 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X \
    --threads 15 --legacy
```

## Evaluate variant calls

Finally, we will evaluate our calls with RTG Tools `vcfeval`. This command requires the reference sequence to be preprocessed:

```
~/tools/rtgtools/rtg format -o ~/data/reference/hs37d5_sdf ~/data/reference/hs37d5.fa
```

Then run vcfeval:

```
rtg vcfeval -t ~/data/reference/hs37d5_sdf \
            -b ~/data/vcf/giab/HG001_GRCh37_truth.vcf.gz \
            --evaluation-regions ~/data/vcf/giab/HG001_GRCh37_hiconf.bed \
            -c ~/data/vcf/NA12878.platinum.b37.octopus.legacy.vcf.gz \
            -o ~/benchmarks/NA12878.platinum.b37.octopus.eval \
            --ref-overlap -f QUAL
```

We see the following results:

```
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
   19.470            3678469        3699288       7611      12392     0.9979       0.9966     0.9973
     None            3679800        3700790       9280      11061     0.9975       0.9970     0.9973
```

## Applying the prototype random forest CSR model

To apply the prototype random forest model, we need to jump through a few more hoop's. First, we need to produce an annotated VCF file to give to the random forest:

```
octopus -R ~/data/reference/hs37d5.fa \
    -I ~/data/bam/NA12878.platinum.b37.bam \
    --filter-vcf ~/data/vcf/NA12878.platinum.b37.octopus.vcf.gz \
    -o ~/data/vcf/NA12878.platinum.b37.octopus.CSR.annotated.vcf.gz \
    -T 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X \
    --csr-train AF SB GQ DP QD FRF CRF GC URF MQ MQ0 MQD
    --threads 15 --legacy
```

Next we create a VCF file filtered with the random forest. The random forest script and model are both in `octopus/prototypes`:

```
./random-forst-csr.py ~/data/vcf/NA12878.platinum.b37.octopus.CSR.annotated.legacy.vcf.gz \
    ~/data/vcf/NA12878.platinum.b37.octopus.random-forest-CSR.legacy.vcf.gz \
    random_forest_n100_d15_model.pkl --drop_info
tabix ~/data/vcf/NA12878.platinum.b37.octopus.random-forest-CSR.legacy.vcf.gz
```

Finally, re-run vcfeval:

```
rtg vcfeval -t ~/data/reference/hs37d5_sdf \
            -b ~/data/vcf/giab/HG001_GRCh37_truth.vcf.gz \
            --evaluation-regions ~/data/vcf/giab/HG001_GRCh37_hiconf.bed \
            -c ~/data/vcf/NA12878.platinum.b37.octopus.random-forest-CSR.legacy.vcf.gz \
            -o ~/benchmarks/NA12878.platinum.b37.octopus.random-forest-CSR.eval \
            --ref-overlap -f QUAL
```

We see the following results:

```
Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------
    2.010            3681191        3702271       3906       9670     0.9989       0.9974     0.9982
     None            3681191        3702271       3906       9670     0.9989       0.9974     0.9982
```
