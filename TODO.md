# TODO

## Features

- Variant filtering.
- Population model/caller.
- Reference callings.
- Reporting realignments.

## Calling accuracy improvements

- The HMM can go over its band limit and give crazy alignments - we need to catch these cases and fallback to a slower routine.
- The MaskOverlappedSegment read transformation should mask segments equally rather than masking just one.
- IndelErrorModel/SnvErrorModel should be polymorphic and dynamically selected at runtime.
- Phaser should allow conditional (on called genotype phasing), and also read supported phasing.
- Model filtering needs to be more selective - filters some TP when haplotypes get long.

## Runtime performance improvements

- Extend pair HMM for AVX.
- Implement SSE log_exp_calculation.
- Allow Callers to parallelise algorithms if in multithreaded mode.
- Model posteriors calculation is very slow, even for individual caller.
- Variational Bayes model needs rewriting as current implementation is just a prototype.
- In multithreaded mode, if we have too many variants buffered, we should write them to another temporary file.

## Cosmetic

- VcfRecordFactory is horrible and needs refactoring. The entire design are the Call family needs looking at.
- Caller::call needs refactoring into smaller methods.

## To consider

- Remove ReferenceGenome dependency from Haplotype.
- Make Genotype::operator[] version non const so can modify in place.
- Implement iterators for ReadReader.

## Testing

- In dire need of proper unit testing!
- Add regression testing
