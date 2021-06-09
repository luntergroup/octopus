---
id: memory
title: Memory Use
---

#### File based factors

Octopus has been designed to minimise disk accesses by buffering frequently accessed resources in main memory. In other words, rather than making lots of short disk reads and using little main memory, octopus prefers to making few large disk reads, and use more main memory. This is beneficial for a number of reasons, especially in a cluster setting. There are essentially two resources which octopus actively buffers: read data and reference sequence. Both these can be controlled by the user:

* `--max-reference-cache-footprint` (`-X`) controls the reference buffer size.
* `--target-read-buffer-footprint` (`-B`) controls the read buffer size. This is not a hard limit, but for most normal samples octopus will respect this.

Both options are specified in all the standard memory units (e.g. `500Mb`, `2G`, `4GB`, etc). Note both buffer sizes are shared between threads, so `--target-read-buffer-footprint=2Gb` using two threads would mean 1GB per thread. This is important is it determines the size of thread 'jobs', and can therefore have a significant impact on throughput and overall runtime.

#### Other factors

The other factors to consider when optimising memory usage are:

* Multithreading: More threads means more memory overhead. A good starting point is to 'budget' an extra 100MB per additional thread, however, this will also depend on your data and the factors below, so some trial and error may be required.
* Calling model: Simpler models use less memory.
* Model setup: The parameterisation of the calling model can have a large influence on short term memory usage. For example setting `--max-joint-genotypes` (in the trio or population model) to a very large memory may mean requests for large memory blocks.