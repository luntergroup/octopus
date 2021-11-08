---
id: forest
title: Random Forest
---

Octopus provides a powerful way to classify variant calls with random forests using the [Ranger](https://github.com/imbs-hl/ranger) random forest library. Pre-trained random forests are available on [Google Cloud](https://console.cloud.google.com/storage/browser/luntergroup/octopus/forests/?project=parabolic-eon-208710). Currently there are forests for germline and somatic variants. You can easily obtain the forests using the Python install script:

```shell
$ ./scripts/install.py --forests
```

will download the forests into `octopus/resources/forests/`. The provided forests are suitable for calling typical diploid germline and cancer data. They *may* work well for other types of samples but this is untested. 

## Using random forest filtering

### Germline variant random forest filtering

To filter germline variants using the germline random forest, just specify the path to the forest in the `--forest-file` option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam \
    --forest resources/forests/germline.forest
```

All calls will be annotated with the `FORMAT` field `RFGQ`, which is the phred-scaled probability of the `GT` call being correct according to the forest model. Each record will also be annotated with an `INFO` field, `RFGQ_ALL`, which is the phred-scalled probability that all `GT` calls in the record are correct, this probability is used to filter calls (see [`--min-forest-quality`](cli.md#--min-forest-quality)) with the `RF` filter.

### Somatic variant random forest filtering

When calling germline and somatic variants, you need to provide both random forests:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam -N NORMAL \
    --forest resources/forests/germline.forest \
    --somatic-forest resources/forests/somatic.forest
```

However, if you are calling somatic variants only (i.e. using the `--somatics-only` flag), then you just need to provide the somatic forest:

```shell
$ octopus -R hs37d5.fa -I normal.bam tumour.bam -N NORMAL \
    --somatics-only \
    --somatic-forest resources/forests/somatic.forest
```

### Re-filtering an Octopus VCF with random forests

Random forest filtering can be applied to a VCF file produced by Octopus without recalling with the `--filter-vcf` option:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam \
    --filter-vcf unfiltered.vcf.gz \
    --forest resources/forests/germline.forest \
    -o filtered.vcf.gz
```

Note this will overwrite any existing filtering information. Note, the VCF file cannot be a `--legacy` VCF file, and should contain all variant calls that Octopus would make by default (i.e. do not use a VCF produced with the `--somatics-only` option).

## Training random forests

The supplied random forest should perform well for the majority of users, however, it is possible that performance can be improved by training a new forest model on your own data.

In addition to Octopus requirements, training new forest models requires the following:

1. A truth set of variants and high-confidence regions (e.g., GIAB or SynDip).
2. [Snakemake](https://snakemake.readthedocs.io/en/stable/).
3. [RTG Tools](https://www.realtimegenomics.com/products/rtg-tools).

Forest models can - and ideally should - be trained on multiple examples runs (i.e., a run of Octopus). Each example can itself be calls for a single sample, or joint calls for multiple samples. For joint calls, each sample is used to generate training data independently, and a subset of samples can be selected.

Training new forests is done using the bundled Snakemake script [forest.smk](https://github.com/luntergroup/octopus/blob/develop/scripts/forest.smk). Basic usage is like:

```shell
$ snakemake \
    --snakefile forest.smk \
    --configfile config.yaml \
    --cores 20
```

You can of course use any [Snakemake option](https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html). Note that the RTG Tools binary `rtg` must be in your path.

The configuration file (YAML or JSON format) provided to Snakemake specifies the training data and setup. The format of the config file is:

* `truths` is a dictionary of truth sets, the key being some unique label (e.g. `SAMPLE1.truth`). The label does **not** need to correspond to the sample name. Each truth set contains a sub dictionary with `vcf` and `bed` value pairs corresponding to the truth variants and high confidence regions, respectively.
* examples is a list of examples to use. Each example is dictionary with the following fields:
    - `name`: the filename prefix to use for this example. If this is not specified then the name is automatically set by concatenating input BAM names. However, this can result in long filenames that may cause OS errors if many inputs are used.
    - `reference`: the reference fasta file. [**required**]
    - 'reads': A list of BAM files to use. The list may be omitted if a single BAM is used. [default none]
    - `truth`: Either the label of the truth set for single sample calling, or a dictionary of sample-truth pairs, where the key is the sample name (in the BAM file), and the value is one of the keys in `truths`. [**required**]
    - `regions`: A bed file to use for calling (i.e. sets option `--regions-file`). [default none]
    - `options`: Additional options to pass to the octopus command.
    - `tp_fraction`: The fraction of true positive calls to use for training. [default 1]
    - `fp_fraction`: The fraction of false positive calls to use for training. [default 1]
    - `threads`: How many threads to use for analysis in this run.
* `training` is a dictionary of forest training setting:
    - `training_fraction`: the fraction of calls to use for training.
    - `hyperparameters`:the random forest hyperparameters to use for training:
        * `trees`: the number of trees to use (ranger option `--ntrees`).
        * `min_node_size`: minimum size of a node in the tree (ranger option `--targetpartitionsize`).
        * `maxdepth`: maximum depth of each tree (ranger option `--maxdepth`).

A minimal example is:

```yaml
examples:
    -
        reference: /path/to/reference.fa
        reads: /path/to/reads.bam
        truth: GIAB//GRCh38//HG001
```

Here, `truth` is set to a special label `GIAB//GRCh38//HG001`. The script accepts special labels like this for any of the GIAB truth sets - it will download all the required truth data. A more detailed example is: 

```yaml
truths:
    SAMPLE1.truth:
        vcf: /path/to/truth/variants/sample1.vcf.gz
        bed: /path/to/truth/variants/sample1.bed
    SAMPLE2.truth:
        vcf: /path/to/truth/variants/sample2.vcf.gz
        bed: /path/to/truth/variants/sample2.bed
examples:
    -
        reference: /path/to/reference.fa
        reads: /path/to/reads1.bam
        truth: SAMPLE1.truth
        tp_fraction: 0.5
        fp_fraction: 1
        threads: 10
    -
        name: joint
        reference: /path/to/reference.fa
        reads:
            - /path/to/reads1.bam
            - /path/to/reads2.bam
        regions: /path/to/calling/regions.bed
        truth:
            SAMPLE1: SAMPLE1.truth
            SAMPLE2: SAMPLE2.truth
        options: --config /path/to/octopus/octopus.config
        threads: 20
training:
    training_fraction: 0.25
    hyperparameters:
        -
            trees: 200
            min_node_size: 20
```

The default behaviour of the script is to use all calls for training, however if the option `--kind=somatic` is used then only variants called `SOMATIC` are used for training. This can be used to generate somatic forests.