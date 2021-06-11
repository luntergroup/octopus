---
id: targeted
title: Targeted Calling
---

By default, octopus will call all contigs specified in the reference index. However, calling can be restricted to a subset of regions by providing a list of regions to call or skip. Note all input regions are assumed to be zero-indexed. If you're using one-indexed regions then add the `--one-based-indexing` option (applied to **all** input regions).

#### `--regions` (`-T`) and `--skip-regions` (`-K`)

Provide a list of regions directly to the command line to call or not call. The format is `<chr>[:start][-[end]]`. So the following are valid:

  * `chr1`: all of `chr1`.
  * `chr2:10,000,000`: the single position `10000000` in `chr2`.
  * `chr3:5,000,000-`: everything from `chr3:5,000,000` onwards.
  * `chr4:100,000,000-200,000,000`: everything between `chr4:100,000,000` and `chr4:200,000,000`. The interval is half open so position `chr4:200,000,000` is **not** included.

You can provide multiple regions with this option, for example:

```shell
$ octopus -R hs37d5.fa -I NA12878.bam \
    -T 1 2:30,000,000- 3:10,000,000-20,000,000
```

Conversely the `--skip-regions` is for providing a list of regions **not** to call. The format is exactly the same as for `--regions`, so:

```console
$ octopus -R hs37d5.fa -I NA12878.bam -K Y
```

Will call all contigs in the reference index other than `Y`. You can provide both `--regions` and `--skip-regions` together, in which case the complement of `--region` and `--skip-regions` will be used.

:::tip

The `--region` option accepts a special ranged argument in the form `<lhs> to <rhs>` where `lhs` and `rhs` are contig names or positions. The reference index is used to expand the range. For example `chr1 to chr3` is expanded to `chr1` `chr2` `chr3`. You can also specify positional start or end points, e.g., `1:10,000,000 to X`. Only one region range can be specified. If it so happens that one of your reference contigs is called `to` then you cannot use this feature!

:::

#### `--regions-file` (`-t`) and `--skip-regions-file` (`-k`)

These commands accept a file containing a line separated list of regions to call or skip. The regions can either be in the same format as `--regions`, or in BED format.

```shell
$ octopus -R hs37d5.fa -I NA12878.bam -t autosomes.txt
$ octopus -R hs37d5.fa -I NA12878.bam -k avoid.bed
```
