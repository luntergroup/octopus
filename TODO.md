# TODO

## Problems

- [] Cancer genotype filtering is not good.
- [] Model filtering needs to be more selective - filters some TP when haplotypes get long.
- [] The MaskOverlappedSegment read transformation should mask segments equally rather than masking just one.
- [] The HMM can go over its band limit and give crazy alignments - we need to catch these cases and fallback to a slower routine.

## Improvements

- [] Improve variant filtering.
- [] Improve IndelErrorModel/SnvErrorModel and make them polymorphic.
- [] Allow Callers to parrallise algorithms - will need to pass policies.
- [] ReadReader should be able to use iterators.

# Refactoring

- [] VcfRecordFactory is pretty horrible
- [] Caller::call needs refactoring into smaller methods

## Bottlenecks

- [] SSE log_exp_calculation.
- [] Assembler.
- [] Phaser.
- [] Variational Bayes models

## To consider

- [] Remove ReferenceGenome dependency from Haplotype
- [] Make Genotype::operator[] version non const so can modify in place

## Testing

- [] In dire need of proper unit testing!
- [] Add regression testing
