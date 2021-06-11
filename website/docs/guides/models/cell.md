---
id: cell
title: Cell
---

The `cell` calling model is used to call germline and somatic variants in single cell and minibatch cell sequencing data. The model attempts to infer local phylogenies for the cells and accounts for allelic biases and dropout often observed in single cell sequencing data.

## Usage: basic

If all of the samples are single cells and none are control cells:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    -o cells.vcf
```

## Usage: with controls

If the experiment includes control cells (e.g. for tumour-normals) then provide the control cell sample names (`--normal-samples`; `-N`):

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --normal-samples CONTROL1 CONTROL2 ... CONTROLM \
    -o cells.vcf
```

All normal cells are assumed to originate from the root (i.e. founder) node of the phylogeny relating cells, and are therefore assumed to all have the same genotype.

## Usage: with minibatchs

If any of the samples are derived from minibatches of cells then specify high dropout concentrations (`--sample-dropout-concentration`) for these samples:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --sample-dropout-concentration MINIBATCH1=100 MINIBATCH2=100 .. MINIBATCHM=100 \
    -o cells.vcf
```

The argument for each minibatch sample may reflect the number of cells contained in the minibatch; the larger the number of cells, the larger the argument value.

The usual use for minibatch samples is for better controls, in which case the minibatches will be normal samples:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --normal-samples MINIBATCH1 MINIBATCH2 ... MINIBATCHM \
    --sample-dropout-concentration MINIBATCH1=100 MINIBATCH2=100 .. MINIBATCHM=100 \
    -o cells.vcf
```

#### VCF output

There are several annotations included in the VCF output:

| Name        | INFO/FORMAT           | Description  |
| ------------- |:-------------:| :-----|
| SOMATIC      | INFO | Indicates that a somatic mutation was inferred (i.e. the phylogeny contains more than one node). |
| PY      | INFO | The MAP phylogeny inferred for the variant loci. This annotation is only added for SOMATIC calls. |
| PPP      | INFO | Posterior probability (Phred) for the MAP phylogeny. |
| PSPP      | INFO | Posterior probabilities (Phred) that the local phylogeny contains `0`,`1`,... nodes |
| PNAP      | FORMAT | Posterior probabilities (Phred) that this sample is assigned to node ID `0`,`1`,.. in the MAP phylogeny (`PY`). |

##### `PY` notation

The phylogeny is serialised using the following algorithm:

```
def serialise(result, node):
    if (node != NULL):
        result += '(' + str(node.id)
        for child in node:
            serialise(result, child)
        result += ')'
```

The algorithm is called with the root node of the phylogeny `serialise("", ROOT)`. Examples:

```
     0
   /   \
  1     2

(0(1)(2))
```

```
     0
       \
        1
         \
          2

(0(1(2)))
```

#### CNV calling

The model can try to identify local copy changes (i.e. deletions or gains of haplotypes). This will result in some samples having called genotypes with different ploidies to the default ploidy. The maximum number of gains and losses is specified with the `--max-copy-gain` and `--max-copy-loss` options, respectively. For example, to identify up to one copy gain or loss:

```shell
$ octopus -C cell \
    -R ref.fa \
    -I cell1.bam cell2.bam ... cellN.bam \
    --max-copy-gain 1 --max-copy-loss 1 \
    -o cells.vcf
```
**Warning** calling copy gains is currently computationally very expensive.

#### Performance considerations

A critical parameter for this calling model is the maximum size of the phylogeny (`--max-clones`). Copy loss and gain calling are also computationally expensive.

It is recommended to allow automatic thread usage with this calling model (use `--threads` option without an argument).