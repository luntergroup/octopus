---
id: errorModels
title: Error Models
---

Octopus accounts for SNV and indel sequencing errors with a context aware error model. The parameterisation of this model is conditional on the library preparations and sequencing technology used, and can have consequences on calling accuracy, particular for indel errors in tandem repeat regions. Octopus comes packaged with parameter sets for several common library preparation and sequencing combinations, and also allows custom sequence error models to be used.

Built-in error models are selected using the `--sequence-error-model` option, which accepts inputs of the form `[library preparation]<.sequencer>`. library preparation is selected from: `PCR`, `PCR-FREE`, or `10X`. sequencer is selected from: `HISEQ-2000`, `HISEQ-2500`, `HISEQ-4000`, `X10`, `NOVASEQ`, `BGISEQ-5000`. For example, `PCR.NOVASEQ` would select the sequence error model parametrised for a `PCR` library preparation and a `NOVASEQ` sequencer. If no sequencer is provided then the default is used (see `octopus --help`).

Custom error models can be used by providing a path to a valid Octopus error model file. These can be produced using the [`profiler.py`](https://github.com/luntergroup/octopus/blob/develop/scripts/profiler.py) Python script in the scripts top level directory. The script creates error model files given the output of the `--data-profile` command line option.
