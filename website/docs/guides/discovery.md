---
id: discovery
title: Variant Discovery
---

The maximum allowed candidate variant size is `--max-variant-size`.

## Pileups

The pileup candidate generator discovers variants directly from input (CIGAR) alignments. Mismatch and indel positions are aggregated for all filtered input reads, and unique candidates are accepted based on several conditions.

The pileup candidate generator can be disabled with the `--disable-pileup-candidate-generator` option.

## Repeat expander

Similar t o

The repeat expander candidate generator can be disabled with the `--disable-repeat-candidate-generator` option.

## Local reassembly



* `--assemble-all`
* Activate regions are tiled into windows up to `--max-assembly-region-size` in length and `--max-assembly-region-overlap` overlap.
* Mismatched read bases with base quality `--assembler-mask-base-quality` are converted to reference bases before threading into the assembly graph. 
  
The assembly graph is then refined to remove spurious bubbles. Bubbles where all wedges have less than `--min-kmer-prune` are removed from the graph. Cycles are removed, unless `--allow-cycles` is specified. 

:::tip

The maximum assembly region size (`--max-assembly-region-size`) implicitly limits the maximum discoverable deletion and complex variant size. However, using larger region sizes increases the chance of getting a reference cycle in the assembly graph, forcing larger k-mer sizes to be used and decreasing sensitivity. 

:::

The assembly candidate generator can be disabled with the `--disable-assembly-candidate-generator` option.

## Input VCF

A set of previously called variants in VCF can also be included as candidates with the `--source-candidates` and `--source-candidates-file` options. By default, only records with `FILTER` set to `PASS` or `.` are included. To use all variants the `--use-filtered-source-candidates` should be specified.

:::important

Supplying input candidates does not guarantee that they will be reported in the final VCF - they are treated just like candidates discovered from the other generators. Octopus does not currently provide complete 'regenotyping' functionality.

:::

:::tip

The default behaviour is to supplement candidates discovered from input reads with input variants. However, it is possible to only use input variants as candidates (e.g. with the `--disable-denovo-variant-discovery` option), which can be useful for certain problems.

:::


