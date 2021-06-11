---
id: preprocessing
title: Read Preprocessing
---

The basic idea of read pre-processing is to remove or modify reads containing artifact sequences that are likely to mislead variant calling or increase runtime by introducing spurious candidates.

## Deduplication

Read duplicates can arise during sequencing in several ways, and are library-prep and technology dependent:

![Docusaurus](/img/guides/duplicates.png)

See [here](https://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/) and [here](http://core-genomics.blogspot.com/2016/05/increased-read-duplication-on-patterned.html) for more detailed discussions on how duplicates arise. 

Duplicates can be problematic for variant calling as they can introduce systematic error (e.g. copying errors during PCR). Removing them is usually recommended for WGS libraries, but this remains somewhat [controversial](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1097-3).




