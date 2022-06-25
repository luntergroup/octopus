---
id: haplotypes
title: Haplotype Proposal
---

The set [discovered candidate variants](guides/../discovery.md) are used to assemble candidate haplotypes in an semi-exhaustive manner. Unlike other haplotype-based variant calling tools such as GATK HaplotypeCaller, Octopus does not require direct read evidence of a haplotype for it to be proposed. The main benifit of this approach is that transative physical read phasing can be exploited to cofidently call haplotypes longer than a single template length. The main downside is the additional runtimes complexity associated with exploring all 2<sup>N</sup> possible haplotypes given N/2 candidate non-overlapping variants. 

## Growth

During the growth phase, novel alleles as appended to the branches of the haplotype tree. 

## Lagging

## Backtracking