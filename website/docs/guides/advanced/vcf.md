---
id: vcf
title: VCF Format
---

Octopus outputs variants in the [VCF 4.3 format](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).

#### Phasing

Octopus always reports phased genotypes (`GT` separated with `|` rather than `/`). The extent of phasing is provided in the `PS` `FORMAT` field, which refers to the `POS` of a previous record, so:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	.	A	C	.	.	.	GT:PS	1|0:100
1	200	.	G	T	.	.	.	GT:PS	0|1:100
```

indicates the variants at positions `100` and `200` are phased. While

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	.	A	C	.	.	.	GT:PS	1|0:100
1	200	.	G	T	.	.	.	GT:PS	0|1:200
```

indicates the variants are unphased. Note that phase sets may not be contiguous.

#### Overlapping sites

To represent complex multi-allelic loci, Octopus prefers a decompose alleles into multiple records and use the `*` allele to resoolve conflicts. In particular, Octopus always splits variants that have unique `REF` alleles into multiple VCF records. For example, two overlapping deletions are represented like:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	102738191	.	ATTATTTAT	A,*	.	.	.	GT	1|2
1	102738191	.	ATTATTTATTTAT	A,*	.	.	.	GT	2|1
```

in contrast to how such a site would be represented by many other tools, either:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTAT	A	.	.	.	GT	1/0
1	102738191	.	ATTATTTATTTAT	A	.	.	.	GT	1/0
```

which is inconsistent as the reference is deduced in each record, or:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTATTTAT	ATTAT,A	.	.	.	GT	1/2
```

which is at least consistent, but rapidly becomes unmanageable as the length and number of overlapping variants increases.
