# Todo

## Features

- [] Finish the assembler.
- [] Multithread everything.

## Performance

- [] Make HMM Verterbi routine faster.
- [] Haplotype & HaplotypeTree currently store Allele's, which contain GenomicRegions, it is redundant
   to store the contig name in each allele as haplotypes must be on a single contig. We should make
   Allele a template class that can use either ContigRegion or GenomicRegion, and use Allele<ContigRegion>
   with Haplotype to save some space.
- [] Genotype currently stores local Haplotype copies. It may be more efficient to store pointers to
   Haplotypes if they will always outlive the genotypes. Investigate this.
- [] We may be able to speed things up by eliminating some haplotypes with an initial fast liklihood
   calculation before entering the EM procedure.
- [] SSE log_exp_calculation.
