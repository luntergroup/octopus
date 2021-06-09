---
id: configs
title: Config files
---

The `--config` command line option is a handy way configure Octopus. The argument is a text file containing settings for any subset of options (other than `config` itself). Each line of the file contains a parameter setting in the format

```shell
parameter = argument
```

Comment lines are allowed in the config file are proceeded with `#`.

It is perfectly fine to specify a config file and explicit command line options, e.g.:

```shell
$ octopus --config my-config.config -R reference.fa -I reads.bam -o octopus.vcf
```

Explicit command line options that are in the config file are ignored.

The [configs](https://github.com/luntergroup/octopus/tree/develop/configs) directory in the main source tree will be used to store useful config files.