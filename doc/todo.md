# Todo

## Problems

- [] Cancer genotype filtering is not good.
- [] ReadIndelErrorModel needs improving.

## Features

- [] Add variant filtering.
- [] Make ReadIndelErrorModel polymorphic.
- [] Multithread everything.

## Performance

- [] SSE log_exp_calculation.
- [] Improve assembler implementation.
- [] ReadReader and VcfReader/Writer should use iterators.
- [] VB models are too slow.
- [] Phasing is too slow.
- [] Assembler is too slow.

## Misc

- [] Remove ReferenceGenome dependency from Haplotype
- [] Make Genotype::operator[] version non const so can modify in place
