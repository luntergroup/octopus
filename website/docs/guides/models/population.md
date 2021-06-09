---
id: population
title: population
---

The `population` calling model is for jointly calling germline variation in small cohorts of samples with known ploidy but unknown pedigree structure. 

## Usage

Multiple samples from the same population, without pedigree information, can be called jointly:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam
```

Joint calling samples may increase calling power, especially for low coverage sequencing.

## Setting ploidies

Octopus can currently joint call samples that have the same ploidies for all contigs. You can set ploidies with the `--organism-ploidy` (`-P`):

```shell
$ octopus -R hs37d5.fa -I haploid1.bam haploid2.bam haploid3.bam -P 1
```

and `--contig-ploidies` (`-p`) options

```shell
$ octopus -R hs37d5.fa -I NA12878.bam NA12891.bam NA12892.bam -p X=1
```

Octopus sets the organism ploidy to 2 by default and also sets chromosome Y to ploidy 1. X is left as two, so if your samples are all male, you could set `-p X=1`.

## Performance considerations

For joint calling, the single most important parameter which will determine performance is [`--max-genotype-combinations`](https://github.com/luntergroup/octopus/wiki/Command-line-reference#option---max-genotype-combinations). This parameter determines the maximum number of joint genotype combinations that octopus can consider. The number of possibl genotype combinations at a given site is exponential in the number of samples, so this value will usually be significantly lower than the number of possible joint genotypes. Octopus will select the 'best' `--max-genotype-combinations` genotype combinations using approximation. Increasing this value may improve accuracy, but result in longer runtimes and memory use.

## Calling large cohorts, and the n + 1 problem

Joint calling many samples (> 10), especially high coverage ones, is computationally expensive. This is because true joint calling necessarily requires considering all possible genotypes for each sample, which grow approximately linearly with the number of samples (due to noise). The number of possible genotype combinations is also polynomial in the number of possible genotypes.

There is also the so-called 'n + 1 problem' of adding an extra sample to a cohort already jointly called. Ideally this would be done cheaply rather than having to recall all samples again. Other methods have addressed this problem by individually calling all samples, then jointly calling all samples by evaluating previously computed genotype likelihoods under a joint model. However, this approach will never correctly call cases where a true allele was not proposed for a specific individual (since no genotype likelihood will be computed). This method also suffers the same 'representation problems' as typical merge procedures as it does not consider haplotypes.

Although octopus does not specifically address the n + 1 problem, we suggest a procedure for jointly calling many samples:

1. Call each sample individually.
2. Group the samples into small subsets (< 20), joint call these using only the variants called previously.
3. Repeat step 2 until there is only one group left (containing all samples).

Adding another sample is achieved in a similar way:

1. Call the new sample individually.
2. Joint call all samples using the variants called in previous joint calling and the new sample.

The reason these procedures is cheaper than joint calling de-novo is the number of candidate alleles is likely significantly reduced as most noise should be removed during individual calling.

For example, three samples could be jointly called with this procedure as follows:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -o octopus.NA12878.bcf
$ octopus -R hs37d5.fa -I NA12891.bam -o octopus.NA12891.bcf
$ octopus -R hs37d5.fa -I NA12892.bam -o octopus.NA12892.bcf
$ octopus -R hs37d5.fa -I NA12892.bam NA12891.bam NA12892.bam \
 --disable-denovo-variant-discovery \
 -c octopus.NA12878.bcf octopus.NA12891.bcf octopus.NA12892.bcf \
 -o octopus.joint.bcf
```

Note the last step disables candidate generation from raw alignments (`-g`) and local reassembly (`-a`), and only uses provided source variants (`-c`).