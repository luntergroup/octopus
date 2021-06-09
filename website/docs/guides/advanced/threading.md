---
id: threading
title: Multithreading
---

#### Usage

Octopus has built in multithreading capabilities, just add the `--threads` command:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads
```

This will let octopus automatically decide how many threads to use (usually the number of available cores). However, a strict upper limit on the number of threads can also be used:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam --threads 4
```

#### Discussion

Using more than one threads is the simplest way to speed up variant calling. However, there is not always a linear payoff in runtime and the number of threads. Optimising thread throughput is challenging as it is highly data dependent. The best case scenario is to have each thread assigned single tasks of equal complexity so that they all finish simultaneously. Internally, octopus divides the calling region into chunks (genomic regions) with approximately the same number of reads. In particular each chunk has [`--target-read-buffer-memory`](https://github.com/luntergroup/octopus/wiki/Command-line-reference#option---target-read-buffer-memory) / `--threads` worth of reads. Allowing octopus more read memory therefore increases the size of each chunk (assuming constant number of threads). Here are some consequences of larger chunks to consider:

* Less thread management overhead (good).
* Less IO accesses (good).
* Higher chance of unbalanced workload (bad).
* More thread blocking due to IO constraint (bad).

Some experimentation may be required to find the best thread/memory combination for your particular data.