---
id: bamout
title: Realigned BAMs
---

Octopus can generate evidence BAMs which can be helpful to understand why a call has been made. Here's an example:

[[images/evidence_bam.png]]

Two conditions must be satisfied in order to make evidence BAMs:

  * All input BAMs must be single sample.
  * An output file (`--output`) must be specified (i.e. no stdout output).

Evidence BAMs are requested using the `--bamout` option. The argument to `--bamout` changes slightly depending on whether you're calling one or more samples: If you're only calling a single sample then the argument to `--bamout` is a *file* path to output to, e.g:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.vcf --bamout octopus.bam
```

For multiple samples the argument to `--bamout` is a *directory* path, e.g:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.vcf --bamout minibams
```

Evidence BAMs with the same names as the input BAMs will be written to this directory, so this cannot be a directory where any of the input BAMs are.

Octopus adds several useful annotations to realigned reads:

| Name        | Description           |
| ------------- |:-------------|
| `HP`      | A list (`,` separated) of haplotype IDs that the read is inferred to originate from. A haplotype ID, which is zero-indexed, corresponds to column in the `GT` field of the affiliated phased VCF. A haplotype ID indicates that the read was unambiguously assigned to the haplotype, while multiple values indicate that the read could equally well be assigned to any of the listed haplotype.  |
| `MD`      | Reference free alignments. As defined in the [SAM specficiation](https://samtools.github.io/hts-specs/SAMtags.pdf) |
| `md`      | Like MD but alignment is relative to the inferred haplotype rather than the reference (i.e. mismatches are inferred sequencing errors). | 
| `hc`      | The CIGAR alignment to the inferred haplotype. |

Coloring alignments by `HP` is a useful way to visualise realigned evidence BAMs:

[[images/tagged_evidence_bam.png]]

In this example we see reads realigned to a HLA-C region. Since all the variants in this region are contained in the same phase set (calls were produced using `--lagging-level AGGRESSIVE`), we can easily see that all heterozygous variant alleles occur on the same haplotype, other than a single bi-allelic allele. The read coloured green was not assigned uniquely to either called haplotype, and is likely miss-mapped. 

By default, only reads supporting regions containing called variation are realigned. However, Octopus can also copy reads overlapping regions where no variation was called using the `--bamout-type FULL` command. Only primary reads are used for BAM realignment.