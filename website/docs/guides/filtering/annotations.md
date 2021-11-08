---
id: annotations
title: Introduction
---

Octopus calls variants and genotypes using a Bayesian model. Like any model, there are some aspects of real error that are not fully captured in the model which can lead to false calls. We therefore recommended filtering calls to improve precision.

There are currently two approaches available in Octopus for filtering variants:

1. [Hard coded thresholds](guides/../thresholds.md)
2. [Random forests](guides/../forest.md)

Both methods use [annotations](#annotations) computed by Octopus.

The random forest approach is preferred when sufficient training data is available (e.g. typical germline and somatic calling). Hard coded thresholds are appropriate for other cases (e.g. UMI low frequency calling).

## Annotations

Octopus provides various annotations for filtering variants. To list available annotations, use the command 

```shell
$ octopus --help --annotations
```

For example, the current annotations are:

```shell
Name	Kind	Number	Type	Description
AC	FORMAT	A	Integer	"Number of non-reference alleles called for each sample"	
AD	FORMAT	R	Integer	"Empirical allele depth"	
ADP	FORMAT	1	Integer	"Number of reads overlapping the position that could be assigned to an allele"	
AF	FORMAT	R	Float	"Empirical allele frequency (AD / ADP)"	
AFB	FORMAT	R	Float	"Absolute difference between empirical allele frequency (AF) and expected allele frequency given genotype"	
AMQ	FORMAT	R	Integer	"Median mapping quality of reads assigned to each allele"	
ARF	FORMAT	1	Float	"Fraction of reads overlapping the call that cannot be assigned to a unique haplotype"	
BMC	FORMAT	1	Integer	"Number of base mismatches at variant position in reads supporting variant haplotype"	
BMF	FORMAT	1	Float	"Fraction of base mismatches at variant position"	
BMQ	FORMAT	1	Integer	"Median quality of base mismatches in reads assigned to a unique allele"	
BQ	FORMAT	R	Integer	"Median base quality of reads supporting each allele"	
CC	INFO	1	Float	"PP divided by QUAL"	
CRF	INFO	1	Float	"Fraction of clipped reads covering the call"	
DAD	FORMAT	R	Integer	"Number of realigned reads supporting ALT alleles identified as duplicates"	
DAF	FORMAT	R	Float	"Fraction of realigned reads supporting ALT alleles identified as duplicates"	
DC	INFO	1	Float	"Number of reads supporting a de novo haplotype in the normal"	
DENOVO	FORMAT	1	Integer	"DENOVO status of each sample"	
DP	FORMAT	1	Integer	"Number of read overlapping the call"	
DPC	FORMAT	1	Float	"Concordance of allele support from duplicate reads"	
ER	FORMAT	R	Float	"Error rate in supporting reads"	
ERS	FORMAT	R	Float	"Error rate standard deviation in supporting reads"	
FRF	FORMAT	1	Float	"Fraction of reads filtered for calling"	
GC	INFO	1	Float	"GC bias of the reference in a window centred on the call"	
GQ	INFO	A	Integer	"Number of non-reference alleles called for each sample"	
GQD	FORMAT	1	Float	"GQ divided by DP"	
ITV	INFO	A	Flag	"Is the variant a transversion"	
MC	FORMAT	1	Integer	"Number of mismatches at variant position in reads supporting variant haplotype"	
MF	FORMAT	1	Float	"Fraction of reads with mismatches at variant position"	
MHL	FORMAT	.	Integer	"Mean likelihood (Phreds) of reads overlapping the site assigned to each haplotype"	
MP	INFO	1	Float	"Model posterior for this haplotype block"	
MQ	INFO	1	Float	"Mean mapping quality of reads overlapping the call"	
MQ0	INFO	1	Integer	"Number of reads overlapping the call with mapping quality zero"	
MQD	FORMAT	1	Integer	"Maximum pairwise difference in median mapping qualities of reads supporting each haplotype"	
MRC	FORMAT	1	Integer	"Number of reads supporting the call that appear misaligned"	
MRL	FORMAT	1	Integer	"Maximum read length overlapping the site"	
NC	INFO	1	Float	"Fraction of overlapping reads supporting a somatic haplotype in the normal"	
PLN	FORMAT	1	Integer	"Length of the phase block for the call"	
PP	INFO	1	Float	"Call posterior probability"	
PPD	INFO	1	Float	"PP divided by DP"	
QD	INFO	1	Float	"QUAL divided by DP"	
QUAL	INFO	A	Integer	"Number of non-reference alleles called for each sample"	
REB	FORMAT	1	Float	"Probability allele occurs at the end (head or tail) of supporting reads"	
REFCALL	FORMAT	1	Integer	"REFCALL status of each sample"	
RSB	FORMAT	1	Float	"Bias of variant side (head or tail half) in supporting reads"	
RTB	FORMAT	1	Float	"Probability allele occurs in the tail of supporting reads"	
SB	FORMAT	1	Float	"Strand bias of reads based on haplotype support"	
SD	FORMAT	1	Float	"Strand bias of reads overlapping the site; probability mass in tails of Beta distribution"	
SF	FORMAT	1	Float	"Max fraction of reads supporting ALT alleles that are supplementary"	
SHC	FORMAT	1	Integer	"Number of called somatic haplotypes"	
SMQ	FORMAT	1	Integer	"Median mapping quality of reads assigned to called somatic haplotypes"	
SOMATIC	FORMAT	1	Integer	"SOMATIC status of each sample"	
STRL	INFO	1	Integer	"Length of overlapping STR"	
STRP	INFO	1	Integer	"Period of overlapping STR"	
VL	FORMAT	1	Integer	"Maximum length of called alleles"
```