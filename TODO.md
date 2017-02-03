# TODO

Here is a non-exhaustive list of things that need doing. In general, the individual caller is now pretty solid (other than some runtime performance issues). Everything else needs work, in particular the cancer caller which suffers from significant runtime issues.

## Features

- Variant filtering.
- Reference callings.

## Calling accuracy improvements

- The HMM can go over its band limit and give crazy alignments - we should try to catch these cases and fallback to a slower routine.
- The MaskOverlappedSegment read transformation should mask segments equally rather than masking just one.
- Model filtering needs to be more selective - filters some TP when haplotypes get long.

## Runtime performance improvements

- Selectively assemble regions: the assembler is slow and many regions don't need assembly, rather than assembling all regions, try to determine regions likely to contain variation and only assemble these.
- Extend pair HMM for AVX.
- Implement SSE log_exp_calculation.
- Allow Callers to parallelise algorithms if in multithreaded mode.
- Implement heuristic model fit calculation before requesting caller model posteriors, as this can be very slow.
- Variational Bayes model needs rewriting as current implementation is just a prototype. In general the CancerCaller is very slow and needs improving.
- In multithreaded mode, if we have too many variants buffered, we should write them to another temporary file.

## Cosmetic

- VcfRecordFactory is horrible and needs refactoring. The entire design are the Call family needs looking at.
- Caller::call needs refactoring into smaller methods.

## To consider

- Remove ReferenceGenome dependency from Haplotype.
- Make Genotype::operator[] version non const so can modify in place.
- Implement iterators for ReadReader.
- Reporting realignments.

## Testing

- In dire need of proper unit testing!
- Add regression testing
