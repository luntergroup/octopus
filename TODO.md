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
- [] Make use of CoverageTracker in CigarScanner (to calculate required read support).
- [] Phaser should allow conditional (on called genotype phasing), and also read supported phasing.
- [] Tasks should start running immediately rather than waiting for all to be generated.
- [] In multithreaded mode, if we have too many variants buffered, we should write them to another temporary file.

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
