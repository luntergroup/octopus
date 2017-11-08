# TODO

This is a non-exhaustive list of things that need doing.

## Features

- Population joint calling.
- Reference callings.

## Calling accuracy improvements

- The HMM can go over its band limit and give crazy alignments - we should try to catch these cases and fallback to a slower routine.
- Model filtering needs to be more selective - filters some TP when haplotypes get long.
- When aggressively phasing, periodically check that phase remains intact, otherwise the number of good haplotypes increases exponentially.

## Runtime performance improvements

- Extend pair HMM for AVX.
- Implement SSE log_exp_calculation.
- Allow Callers to parallelise algorithms if in multithreaded mode.
- Implement heuristic model fit calculation before requesting caller model posteriors, as this can be very slow.
- Variational Bayes model needs rewriting as current implementation is just a prototype. In general the CancerCaller is very slow and needs improving.
- In multithreaded mode, if we have too many variants buffered, we should write them to another temporary file.

## Refactoring

- VcfRecordFactory is horrible and needs refactoring. The entire design are the Call family needs looking at.

## To consider

- Remove ReferenceGenome dependency from Haplotype.
- Make Genotype::operator[] version non const so can modify in place.
- Implement iterators for ReadReader.
- Reporting realignments.

## Testing

- More unit testing!
- Add regression testing,
