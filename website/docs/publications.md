---
id: publications
title: Publications
---

Here are the current papers associated with Octopus.

## [A unified haplotype-based method for accurate and comprehensive variant calling](https://www.nature.com/articles/s41587-021-00861-3)

This is the main Octopus paper describing the method with benchmarks for germline and somatic variant calling.

#### Abstract

Almost all haplotype-based variant callers were designed specifically for detecting common germline variation in diploid populations, and give suboptimal results in other scenarios. Here we present Octopus, a variant caller that uses a polymorphic Bayesian genotyping model capable of modeling sequencing data from a range of experimental designs within a unified haplotype-aware framework. Octopus combines sequencing reads and prior information to phase-called genotypes of arbitrary ploidy, including those with somatic mutations. We show that Octopus accurately calls germline variants in individuals, including single nucleotide variants, indels and small complex replacements such as microinversions. Using a synthetic tumor data set derived from clean sequencing data from a sample with known germline haplotypes and observed mutations in a large cohort of tumor samples, we show that Octopus is more sensitive to low-frequency somatic variation, yet calls considerably fewer false positives than other methods. Octopus also outputs realigned evidence BAM files to aid validation and interpretation.

#### Bibtex

```tex
@article{octopus,
   author = {Cooke, Daniel P. and Wedge, David C. and Lunter, Gerton},
   title = {A unified haplotype-based method for accurate and comprehensive variant calling},
   journal = {Nature Biotechnology},
   ISSN = {1546-1696},
   DOI = {10.1038/s41587-021-00861-3},
   url = {https://doi.org/10.1038/s41587-021-00861-3
https://www.nature.com/articles/s41587-021-00861-3.pdf},
   year = {2021},
   type = {Journal Article}
}
```

## [Benchmarking small-variant genotyping in polyploids](https://genome.cshlp.org/content/32/2/403)

In this paper, we benchmark Octopus on polyploid samples.

#### Abstract

Genotyping from sequencing is the basis of emerging strategies in the molecular breeding of polyploid plants. However, compared with the situation for diploids, in which genotyping accuracies are confidently determined with comprehensive benchmarks, polyploids have been neglected; there are no benchmarks measuring genotyping error rates for small variants using real sequencing reads. We previously introduced a variant calling method, Octopus, that accurately calls germline variants in diploids and somatic mutations in tumors. Here, we evaluate Octopus and other popular tools on whole-genome tetraploid and hexaploid data sets created using in silico mixtures of diploid Genome in a Bottle (GIAB) samples. We find that genotyping errors are abundant for typical sequencing depths but that Octopus makes 25% fewer errors than other methods on average. We supplement our benchmarks with concordance analysis in real autotriploid banana data sets.

#### Bibtex

```tex
@article{Cooke01022022,
author = {Cooke, Daniel P. and Wedge, David C. and Lunter, Gerton}, 
title = {Benchmarking small-variant genotyping in polyploids},
volume = {32}, 
number = {2}, 
pages = {403-408}, 
year = {2022}, 
doi = {10.1101/gr.275579.121}, 
URL = {http://genome.cshlp.org/content/32/2/403.abstract}, 
eprint = {http://genome.cshlp.org/content/32/2/403.full.pdf+html}, 
journal = {Genome Research} 
}
```

## [Accurate genotyping of single cells with Octopus](https://www.researchsquare.com/article/rs-583831/v1)

#### Absract

We describe an extension to our variant calling tool, Octopus (https://github.com/luntergroup/octopus), for single-cell DNA sequencing data. Octopus jointly genotypes cells from a lineage, accounting for amplification stochasticity and sequencing error with a haplotype-based Bayesian model. Octopus is considerably more accurate at genotyping single cells than existing methods.

#### Bibtex

```tex
@article{RN862,
   author = {Daniel, Cooke P. and Gerton, Lunter and David C., Wedge},
   title = {Accurate genotyping of single cells with Octopus},
   journal = {Research Square},
   ISSN = {2693-5015},
   DOI = {10.21203/rs.3.rs-583831/v1},
   url = {https://doi.org/10.21203/rs.3.rs-583831/v1},
   year = {2021},
   type = {Journal Article}
}

```