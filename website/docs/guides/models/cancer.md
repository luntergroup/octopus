---
id: cancer
title: Cancer
---

The `cancer` calling model is for calling germline variation and somatic mutations in tumours. The model can jointly genotype multiple tumours from the same individual, and make use of a normal sample for improved classification power.

## Usage: tumour-normal pairs

To call germline and somatic mutations in a paired tumour-normal sample, just specify which sample is the normal (`--normal-sample`; `-N`):

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam --normal-sample NORMAL
```

It is also possible to genotype multiple tumours from the same individual jointly:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumourA.bam tumourB.bam -N NORMAL
```

## Usage: tumour only

If a normal sample is not present the cancer calling model must be invoked explicitly:

```shell
$ octopus -R hs37d5.fa -I tumourA.bam tumourB.bam -C cancer
```

Be aware that without a normal sample, somatic mutation classification power is significantly reduced.

## VCF output

By default both germline and somatic variants are called, somatic mutations are tagged with the `SOMATIC` INFO field. The `GT` fields for `SOMATIC` variants and any other variants in the same phase set (`PS`) are augmented with the number of unique somatic haplotypes inferred (somatic haplotypes are identified with the `HSS` tag). For example, in the following VCF:

```
##fileformat=VCFv4.3
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Record includes a somatic mutation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set">
##FORMAT=<ID=HSS,Number=.,Type=Integer,Description="Somatic status for each haplotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NORMAL TUMOUR
1 50  . G  A   . . SOMATIC GT:PS:HSS   0|0:50:0,0    1|0|0:50:1,0,0
1 100 . A  C   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|0|1:100:0,0,1,1
1 150 . G  T   . . .       GT:PS:HSS   1|0:100:0,0   1|0|1|0:100:0,0,1,1
1 200 . A  AC  . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|0:100:0,0,1,1
1 250 . TA T   . . .       GT:PS:HSS   1|1:100:0,0   1|1|1|1:100:0,0,1,1
2 300 . C  A,T . . SOMATIC GT:PS:HSS   0|2:300:0,0   1|0|2:300:1,0,0
2 350 . G  T   . . .       GT:PS:HSS   0|1:300:0,0   0|0|1:300:1,0,0
2 400 . T  C   . . SOMATIC GT:PS:HSS   0|0:300:0,0   1|0|0:300:1,0,0
3 100 . A  C   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|1:100:0,0,1,1
3 150 . C  CC  . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|0:100:0,0,1,1
3 200 . G  T   . . .       GT:PS:HSS   0|1:100:0,0   0|1|1|1:100:0,0,1,1
4 100 . C  T   . . .       GT:PS:HSS   0|1:100:0,0   0|1|1|1:100:0,1,1,0
4 150 . A  G   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|0|1|0:100:0,1,1,0
4 200 . A  G   . . SOMATIC GT:PS:HSS   0|0:100:0,0   0|1|0|0:100:0,1,1,0
```

The first phase set `1:50` includes a simple somatic mutation; it is not phased with any other variants. Downstream of this is another phase block starting at `1:100` that includes 2 germline and 2 somatic variants. The first somatic mutation in this phase set at `1:100` is phased onto the germline haplotype including the reference allele at the germline variant at `1:150`. The second somatic mutation in this phase set is phased with the alternate allele of this germline variant. In the third phase set stating at `2:300` there is a somatic mutation that segregates with a germline variant at the same position, and another somatic mutation which is phased onto the same germline haplotype - the reference. The somatic allele in the multiallelic record is determined by looking at the `HSS` flags (1 indicates somatic), so "C>A" is the somatic mutation in this case. The third phase set starting at `3:100` includes 2 somatic mutations phased onto the same germline haplotype, but the first somatic mutation was inferred to segregate with both the germline allele and somatic allele at `3:150`, suggests linear progression of the C>CC somatic mutation. The last phase set starting at `4:100` includes 2 somatic mutations phased onto the same germline haplotype, but not with each other.

## `QUAL` vs `PP`

For both paired and tumour-only calling, octopus reports two quality scores for each call (both germline and somatic):

  * `QUAL` is the posterior probability the variant is *segregating* in the sample regardless of somatic classification.
  * `PP` (an `INFO` field) is the posterior probability the variant is segregating **and** classified correctly.

The difference between `QUAL` and `PP` indicates the uncertainty in the calls classification; a call may have high `QUAL` but low `PP` if the classification is uncertain (common in tumour-only calling or if the normal coverage is low). `PP` should always be less than `QUAL` in theory.

## `HPC`, `MAP_HF` and `HF_CR`

Octopus infers a probability distribution over haplotype frequencies, including any somatic haplotypes. For each variant in a phase-block with a `SOMATIC` variant, three statistics are reported that relate to haplotype frequency (per haplotype).

  * `HPC` is the Dirichlet posterior pseudo-count. Note that this count includes the prior count so should not be taken to mean the empirical count.
  * `MAP_HF` (`FORMAT`) is the [Maximum a Posteriori](https://en.wikipedia.org/wiki/Maximum_a_posteriori_estimation) haplotype frequency point estimate.
  * `HF_CR` (`FORMAT`) is a [credible interval](https://en.wikipedia.org/wiki/Credible_interval) of the haplotypes frequency. The mass of the credible interval is specified by `--credible-mass`. The credible interval gives you an indication of how certain the MAP_HF estimate is; A very narrow interval means the MAP estimate is very certain, a wide interval means it is uncertain.

## Performance considerations

* The number of genotypes considered by the model `--max-genotype` has a significant impact on overall runtime.
* The parameter `--max-somatic-haplotypes` controls the maximum number of unique segregating somatic haplotypes to be modelled. There must be at-least one somatic haplotype, but adding more can resolve somatic mutations falling on different germline haplotypes or multiple distinct haplotypes due to sub clonal evolution. 

## SOMATICs only

To report only `SOMATIC` calls just add the `--somatics-only` command. This will only work if calls are beings filtered already (i.e. it will have no effect if `-f off`). It is generally not recommended to use this option until you are 100% satisfied with your calls, as call filtering will not work correctly if germline variants have been filtered from the VCF; you will not be able to re-filter your SOMATIC-only VCF.