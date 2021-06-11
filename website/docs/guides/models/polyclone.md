---
id: polyclone
title: Polyclone
---

The `polyclone` calling model is for calling variation in a pooled sample of haploid clones where the number and mixture composition of clones is unknown, as is sometimes the case in bacteria or virus sequencing studies. 

## Usage

If your sample contains an unknown mix of haploid clones (e.g. some bacteria or viral samples), use the `polyclone` calling model:

```shell
$ octopus -R H37Rv.fa -I mycobacterium_tuberculosis.bam -C polyclone
```

This model will automatically detect the number of subclones in your sample (up to the maximum given by `--max-clones`).

## Performance considerations

The most important parameter in this model is `--max-clones` which determines the maximum number of haplotypes that octopus will try to use to fit the data. This is a bit like setting the 'ploidy' for the sample. Higher values of `--max-clones` may lead to longer runtimes and more memory usage - we do not recommend values above 10.