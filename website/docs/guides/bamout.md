---
id: bamout
title: Realigned BAMs
---

Octopus can generate realigned BAMs that provide visual evidence for why a call has been made. Realigned BAMs are particularly helpful for confirming complex variation where the mapper alignments are incorrect, as can be seen in the IGV pileups below 

<iframe width="750" height="600" src="https://igv.org/web/release/2.8.5/embed.html?sessionURL=blob:rVTbitswEP2VRU8teGU7VytQCm3ptlD6sNvSh1LK2B472tiSI8l2LuTfd.x19sKmEMqGxMTM6Mw5ZzSzZw0aK7ViCzbiEZ8yj9mlbm.grAr8DiVatsigsOgxgxkaVAmyxZ7JlE4s83FEBxSl0duXugR18ebq.uNyHPld7C0FM7AOfl5_69Kdq.zC9.2YQwk7raC1PNGlL_OGx0ZDKpV10tUOuTa5n6PSRMC3uO7h.gfPgFClSnHzuqj0lYScbJ2OQaVngBMaP6L1SKCUduDITut3OB8I52uKmruN4_mO0KGQYP8Hunv87U9zBzE7eKzQSU3NYcnShIvIm80ibz4NLrt_wguigKo5A8mKcn7vmdtWXY9IdN230GPapGjY4lIEwTwUYjSdzCeBEOHB27PaFE8Ytm3LU6OrWG96gtYXo1rd5rjTEwHBcuXfbNUnWfG4hcsSSx5D.T4t3nUUhrvxMuF0E1.WMtpUKs_WmRJ6OytPlKKfPJbLtCnBEeB9hUE1.ZarEpXrO1BggZ8Nrn8sDdJdL.gmB3z0YMjkHAPW8rax2Xwz1jJcResjK504XdX23wY8STjXAJxGtRgF1Bqduiw5UeqVDZieY4COZ5GATbJp1Dhah0dWTZLRRT.t_T52rmyx263Kqhq7KhcwC54X4C5.qZhCj4obMBJ6vYOqWadqIHSNGQ3CxRUqWm9PEGjD0dB1w_Fc_nkjqpJYEnKH.jjwvVYkh4cl2kgrY1lIt_1FId3SAIbdai11A3FB5Ia8gXUY9J_HRj5sGHb4c7gD">
</iframe><br /><br />

Evidence BAMs are requested using the `--bamout` option. The argument to `--bamout` changes slightly depending on whether you're calling one or more samples: If you're only calling a single sample then the argument to `--bamout` is a file path to write the BAM to, e.g.:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.vcf --bamout octopus.bam
```

For multiple samples the argument to `--bamout` is a directory path, e.g.:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.vcf --bamout minibams
```

Realigned BAMs with the same names as the input BAMs will be written to this directory, so this cannot be a directory where any of the input BAMs are located.

:::important

Realigned BAMs are only available for single-sample BAMs and when `--output` is specified (i.e. no stdout output).

:::

Octopus adds several useful annotations to realigned reads:

| Name        | Description           |
| ------------- |:-------------|
| `HP`      | A list (`,` separated) of haplotype IDs that the read is inferred to originate from. A haplotype ID, which is zero-indexed, corresponds to column in the `GT` field of the affiliated phased VCF. A haplotype ID indicates that the read was unambiguously assigned to the haplotype, while multiple values indicate that the read could equally well be assigned to any of the listed haplotype. |
| `MD`      | Reference free alignments. As defined in the [SAM specficiation](https://samtools.github.io/hts-specs/SAMtags.pdf) |
| `md`      | Like MD but alignment is relative to the inferred haplotype rather than the reference (i.e. mismatches are inferred sequencing errors). | 
| `hc`      | The CIGAR alignment to the inferred haplotype. |
| `PS`      | The phase set the read was assigned to. |

:::tip

The `HP` tag is useful for colouring and grouping alignments in IGV.

:::

By default, only reads supporting regions containing called variation are realigned. However, Octopus can also copy reads overlapping regions where no variation was called using the `--bamout-type FULL` command. Only primary reads are used for BAM realignment.

:::caution

Reads are assigned and realigned to haplotypes called in the `--output` VCF. This means that read-pairs in different phase sets can appear discordant, and reads that are not completely spanned by a phase set (or overlap multiple phase sets) may have poor alignments. Consider trying to [increase haplotype lengths](phasing.md) if this occurs.

:::