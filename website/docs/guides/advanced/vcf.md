---
id: vcf
title: VCF Format
---

Octopus reports variants in [VCF 4.3 format](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.pdf).

## Haplotypes

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

## Spanning alleles

To represent complex multi-allelic loci, Octopus prefers a decompose alleles into multiple records and use the `*` allele to resolve conflicts. In particular, Octopus always splits variants that have unique `REF` alleles into multiple VCF records. For example, two overlapping deletions are represented like:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	102738191	.	ATTATTTAT	A,*	.	.	.	GT	1|2
1	102738191	.	ATTATTTATTTAT	A,*	.	.	.	GT	2|1
```

In contrast, some other tools may choose to represent these variant with

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTAT	A	.	.	.	GT	1/0
1	102738191	.	ATTATTTATTTAT	A	.	.	.	GT	1/0
```

which is inconsistent as the reference is deduced in each record, or

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	102738191	.	ATTATTTATTTAT	ATTAT,A	.	.	.	GT	1/2
```

which is consistent, but rapidly becomes unmanageable as the length and number of overlapping variants increases.

### Multiple `*` records

In some cases, a site may overlap with multiple distinct alleles. Octopus represents these sites using multiple `*` alleles in the `ALT` field.

Consider the [realignments](../bamout.md) from the [Ashkenazim trio](https://www.coriell.org/1/NIGMS/Collections/NIST-Reference-Materials):

<iframe width="750" height="700" src="https://igv.org/web/release/2.8.5/embed.html?sessionURL=blob:rZRva9swEMa_StGrDVzZcZLZDozBBmsHYy.yjb0YY8jy2VYiS64k2_lDvnvPatJ0bIUwGhKTSKdfnufudHvSg7FCK7IgMU3pnATE1nr4yppWwhfWgCWLkkkLATFQggHFgSz2RBR4oq6mKR5QGIa_bruGqatXN8sP9TQNx73XuFky69j35ecx3LnWLsLQTilr2E4rNljKdROKqqe50awQyjrhOgdUmyqsQGkUEFq48zj_oCVDqlAFbF6Wim.BZL51OmequACONHqieRJTSjvmMJ02HDnvkfOpAE3dxtFqh3QmBbP_gx4fv_1p6lhODgGRmndYHMJrM1tMsiCexkESpdenr2mC_.cM42uM.rknbtuOVULbnS9iQLQpwJDFdRZFySTL4vksmUVZNjkEe9IZ.UTjMAy0MLrN9cZLtOFGDNtVtlt1nZg3fR_e3kRRTHPWvCvk2.hJU5zW_120v8Hrsl8XQ5mUd0NVRuYMxo84wUttGuaQ8wA.WsP0VKoB5XyiJUj4aODuW20AW1piw0Y0fnQ9u8RllXUy1jYxepOUzYPL6TMu_fqlLlet4Mk0Wsluu2tbdga_sMv5JS53XaPrbcpz3oGsuRcze8alX7_Upd7VrQFd5yu.XsflGfzCLt.MLo8il1Bii1_dgMLR9YSO0wsv1Nj2f6bjsuuneC6QPFLPl9nnAFDLcUD2wopcSOG2P3BLD3i1JuPYbHTPconijnFH1ZPIv86WH6cHOfw63AM-">
</iframe><br /><br />

The variants called by Octopus in this region are:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002	HG003	HG004
chr4	19232687	.	A	G	50	PASS	AC=2;AN=6	GT:PS	0|1:19232687	0|0:19232687	0|1:19232687
chr4	19232728	.	ATCTG	A	50	PASS	AC=2;AN=6	GT:PS:PQ	0|1:19232687	0|0:19232687	0|1:19232687
chr4	19232732	.	GTCTGTCTATCTA	G,*,*	50	PASS	AC=2,2,2;AN=6	GT:PS	2|3:19232687	1|2:19232687	1|3:19232687
chr4	19232736	.	G	GTCTATCTA,*	50	PASS	AC=2,2;AN=6	GT:PS	1|0:19232687	2|19232687	2|0:19232687
chr4	19232762	.	C	CTA	550	PASS	AC=2;AN=6	GT:PS	0|1:19232687	0|0:19232687	0|1:19232687
chr4	19232762	.	CTATATA	C	50	PASS	AC=2;AN=6	GT:PS	0|0:19232687	1|0:19232687	1|0:19232687
chr4	19232858	.	A	ACT	50	PASS	AC=4;AN=6	GT:PS	0|1:19232687	1|0:19232687	1|1:19232687
```

Notably, Octopus reports two `*` alleles at the `chr4:19232732` record. This indicates that there are two distinct alleles overlapping this site (at `chr4:19232728` and `chr4:19232736`). 

:::caution

In some rare cases, Octopus may report double `*` records in a single sample, resulting in calls like

```
1	100	.	ACGT	*,*	50	PASS	AC=1,1;AN=2	GT	2|3
```

This generally indicates a homozygous deletion directly upstream of a heterozygous deletion, like:

```
1	99	.	TA	T	50	PASS	AC=2;AN=2	GT:PS	1|1:99
1	100	.	ACGT	*,*	50	PASS	AC=1,1;AN=2	GT:PS	2|3:99
```

:::

