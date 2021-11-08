---
id: cli
title: Command Line Reference
---

## General options

### `--help`

Command `--help` (short `-h`) prints the list of available command line options.

```shell
$ octopus --help
Octopus command line options:

General:
  -h [ --help ]                         Report detailed option information
  --version                             Report detailed version information
  --config arg                          Config file to populate command line 
                                        options
  --options                             Log all command line option values at 
                                        startup
  --debug [=arg(="octopus_debug.log")]  Create log file for debugging
  --trace [=arg(="octopus_trace.log")]  Create very verbose log file for 
                                        debugging
  -w [ --working-directory ] arg        Sets the working directory
  --threads [=arg(=0)]                  Maximum number of threads to be used. 
                                        If no argument is provided unlimited 
                                        threads are assumed
<OUTPUT STRIPPED>
```

### `--version`

Command `--version` prints the version of the current binary, including the git branch and commit. The command also prints some system information about the binary.

```shell
$ octopus --version
octopus version 0.6.3-beta (develop e7205439)
Target: x86_64 Darwin 18.7.0
SIMD extension: AVX2
Compiler: AppleClang 10.0.1.10010046
Boost: 1_70
```

### `--config`

Option `--config` is used to specify [configuration files](guides/advanced/configs.md).

```shell
$ octopus ---config myconfig.txt
```


### `--trace`

Option `--trace` results in a log file containing detailed debug information being created.

```shell
$ octopus -R ref.fa -I reads.bam --trace # default name is "octopus_trace.log"
$ octopus -R ref.fa -I reads.bam --trace =/foo/bar/trace.log
```

**Warning** Trace files can get very large, so only use on small inputs..

### `--working-directory`

Option `--working-directory` (short `-w`) is used to set the working directory of the run. All output and temporary files will be relative to the working directory, unless absolute file paths are provided.

```shell
$ octopus -w ~/vcf -o octopus.vcf
```

The output will be in `~/vcf/octopus.vcf`.

**Notes**

* If no working directory is specified then the current directory is used.

### `--resolve-symlinks`

Command `--resolve-symlinks` forces all symbolic link input paths to be replaced with their resolved targets during program initialisation.

```shell
$ octopus -R ref.fa -I reads.bam --resolve-symlinks
```

### `--threads`

Option `--threads` is used to enable [multithreaded execution](guides/advanced/threading.md). The option accepts a positive integer argument that specifies the maximum number of threads to use:

```shell
$ octopus -R ref.fa -I reads.bam --threads 4 # use maximum of 4 threads
```

If no argument is provided then automatic thread handling is used. This may be greater than the number of system threads.

```shell
$ octopus -R ref.fa -I reads.bam --threads #automatic thread handling
```

### `--max-reference-cache-memory`

Option `--max-reference-cache-memory` (short `-X`) controls the size of the buffer used for reference caching, and is therefore one way to [control memory use](guides/advanced/memory.md). The option accepts a non-negative integer argument in bytes, and an optional unit specifier.

```shell
$ octopus -R ref.fa -I reads.bam -X 0 # disable reference caching
$ octopus -R ref.fa -I reads.bam -X 1g # reference cache is 1 gigabyte
```

**Notes**

* Capitisation of the units is ignored.
* An argument of `0` disables reference caching.
* The buffer size never exceeds the size of the reference.

### `--target-read-buffer-memory`

Option `--target-read-buffer-memory` (short `-B`) controls the size of the memory buffer used for input sequencing reads, and is therefore one way to [control memory use](guides/advanced/memory.md). The option accepts a positive integer argument in bytes, and an optional unit specifier.

```shell
$ octopus -R ref.fa -I reads.bam -B 20G # reference cache is 20 gigabyte
```

**Notes**

* Capitisation of the units is ignored.
* The minimum read buffer size is `50Mb`; arguments less than this are ignored.

### `--target-working-memory`

Option `--target-working-memory` sets the target amount of working memory for computation, and is therefore one way to [control memory use](guides/advanced/memory.md). The option accepts a positive integer argument in bytes, and an optional unit specifier. The option is not strictly enforced, but is sometimes used to decide whether to switch to lower-memory versions of some methods (possibly at the cost of additional runtime).

```shell
$ octopus -R ref.fa -I reads.bam --target-working-memory 40Gb
```

**Notes**

* Capitisation of the units is ignored.

### `--temp-directory-prefix`

Option `--temp-directory-prefix` sets the Octopus working temporary directory prefix.

```shell
$ octopus -R ref.fa -I reads.bam --temp-directory-prefix ~/octopus_tmp
```

**Notes**

* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.
* If a directory with the prefix already exists, then `-N` is appended to the prefix where `N` is the smallest integer such that the resulting path does not exist.

### `--reference`

Option `--reference` (short `-R`) sets the reference FASTA file used for variant calling. This is a required option.

```shell
$ octopus -R ref.fa -I reads.bam
```

**Notes**

* In addition to the FASTA file, a reference FASTA index file (extension `.fai`) is also required. This must be present in the same directory as the FASTA file.
* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.

### `--reads`

Option `--reads` (short `-I`) specifies the input read files (extension `.bam` or `.cram`) used for variant calling. The argument is a list of file paths. This is a required option (unless `--reads-file` is specified).

```shell
$ octopus -R ref.fa -I reads1.bam reads2.bam
```

**Notes**

* Each read file must have a paired index file in the same directory.
* The option can be specified multiple times on the command line (arguments are concatenated).
* Can be used in conjunction with `--reads-file` (arguments are concatenated).
* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.


### `--reads-file`

Option `--reads-file` (short `-i`) specifies a file containing a list of input read files (one per line).

```shell
$ octopus -R ref.fa -i bams.txt
```

### `--regions`

Option `--regions` (short `-T`) specifies a list of genomic [regions to call](guides/advanced/targeted.md).

```shell
$ octopus -R ref.fa -I reads.bam -T chr1 # all of chr1
$ octopus -R ref.fa -I reads.bam -T chr1 chr2 # all of chr1 and chr2
$ octopus -R ref.fa -I reads.bam -T chr1:0-1,000 # first 1,000bp of chr1
$ octopus -R ref.fa -I reads.bam -T chr1:100- # everything after position 100 of chr1
$ octopus -R ref.fa -I reads.bam -T chr1 to chrM # as specified in reference index
```

**Notes**

* If this option is not specified (and neither is `--regions-file`), then all contigs present in the reference genome index are used.
* Commas in the position tokens are ignored.
* The option can be specified multiple times on the command line (arguments are concatenated).

### `--regions-file`

Option `--regions-file` (short `-t`) specifies a file containing a list of genomic regions to call (one per line). The file format can either be plain text (using the same input format as the `--regions` option), or BED format (i.e. tab-separated). GZipped files are also accepted.

```shell
$ octopus -R ref.fa -I reads.bam -t regions.bed
```

### `--skip-regions`

Option `--skip-regions` (short `-K`) specifies a list of genomic regions to ignore during calling.

```shell
$ octopus -R ref.fa -I reads.bam -K X Y MT
```

**Notes**

* The option can be specified multiple times on the command line (arguments are concatenated).

### `--skip-regions-file`

Option `--skip-regions-file` (short `-k`) specifies a file containing a list of genomic regions to skip (one per line). The file format can either be plain text (using the same input format as the `--regions` option), or BED format (i.e. tab-separated). GZipped files are also accepted.

### `--one-based-indexing`

Command `--one-based-indexing` directs the program to read all input regions (i.e. those specified in options `--regions`, `--regions-file`, `--skip-regions`, and `--skip-regions-file`) using 1-based indexing rather than 0-based indexing.

```shell
$ octopus -R ref.fa -I reads.bam -T chr1:100 --one-based-indexing
$ octopus -R ref.fa -I reads.bam -T chr1:99 # equivalent to above
```

### `--ignore-unmapped-contigs`

Command `--ignore-unmapped-contigs` can be used to force execution if there is a mismatch between the input reference genome and the one used to map input reads (as specified in the SAM header). In particular, if there are contigs present in the reference genome but not in the reference used for read mapping.

```shell
$ octopus -R ref.fa -I reads.bam
<INFO> ------------------------------------------------------------------------
<INFO> octopus v0.6.1-beta (develop 9c97253c)
<INFO> Copyright (c) 2015-2019 University of Oxford
<INFO> ------------------------------------------------------------------------
<EROR> A user error has occurred:
<EROR> 
<EROR>     Some or all of the contigs in the reference genome (ref) are not
<EROR>     present in the read files.
<EROR> 
<EROR> To help resolve this error ensure the reference genome used for mapping
<EROR> is the same as the one used for calling (ref) and all input contigs are
<EROR> present in the read headers.
<INFO> ------------------------------------------------------------------------

$ octopus -R ref.fa -I reads.bam --ignore-unmapped-contigs # ignore error above

```

### `--samples`

Option `--samples` (short `-S`) specifies a list of samples to use for variant calling. Each sample must be present in the input read files (using the `SM` SAM tag).

```shell
$ octopus -R ref.fa -I reads.bam -S sample1 sample2
```

**Notes**

* If no samples are specified, all samples present in the input read files are used.
* The option can be specified multiple times on the command line (arguments are concatenated).

### `--samples-file`

Option `--samples-file` (short `-s`) specifies a file containing a list (one per line) of samples to use for variant calling.

```shell
$ octopus -R ref.fa -I reads.bam -s samples.txt
```

### `--output`

Option `--output` (short `-O`) sets the output destination. The file extension is used to determine the output file type (valid types are `.vcf`, `.vcf.gz`, and `.bcf`). If the file type is not recognised then uncompressed VCF is used.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.bcf
```

**Notes**

* If no output is specified, `stdout` is used.
* Like all path arguments in Octopus, the argument is assumed to be relative to the `--working-directory`, unless the path is absolute and exists.

### `--contig-output-order`

Option `--contig-output-order` specifies the order that records will be processed and written to the output. Possible options are: `lexicographicalAscending`, `lexicographicalDescending`, `contigSizeAscending`, `contigSizeDescending`, `asInReferenceIndex`, `asInReferenceIndexReversed`, `unspecified`

```shell
$ octopus -R ref.fa -I reads.bam --contig-output-order asInReferenceIndexReversed
```

### `--sites-only`

Command `--sites-only` removes genotype information from the final output (i.e. drops the VCF `FORMAT` and sample columns).

```shell
$ octopus -R ref.fa -I reads.bam --sites-only
```


### `--bamout`

Option `--bamout` is used to produce [realigned evidence BAMs](guides/bamout.md). The option input is a file prefix where the BAMs should be written.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.vcf --bamout reads.realigned.bam # if 'reads.bam' contains one sample
$ octopus -R ref.fa -I readsA.bam readsB.bam -o calls.vcf --bamout minibams
```

**Notes**

* The `--bamout` command is only allowed if the final output is written to file (i.e. with `--output`).
* The default behaviour of `--bamout` is to realign only those reads supporting variant haplotypes (see also `--full-bamout`).


### `--bamout-type`

Option `--bamout-type` is used to select the type of realigned evidence BAM to produce when using the `--bamout` option. If `FULL` is chosen then all reads in the input BAM are written to the evidence BAM. If `MINI` is chosen (default) then only reads supporting variant sites are written.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.vcf --bamout reads.realigned.bam --bamout-type FULL
```

### `--pedigree`

Option `--pedigree` is used to input a [pedigree file](http://www.helsinki.fi/~tsjuntun/autogscan/pedigreefile.html) (`.ped`). Some calling models may use pedigree information to improve calling accuracy.

```shell
$ octopus -R ref.fa -I mother.bam father.bam child.bam --pedigree family.ped 
```

where `family.ped` might contain, for example

```text
my_trio	FATHER	0	0	1	0
my_trio	MOTHER	0	0	2	0
my_trio	CHILD	FATHER	MOTHER	2	0
``` 

**Notes**

* Specifying this option with a `.ped` file defining a trio relationship between three input samples will automatically activate the `trio` calling model.

### `--fast`

Command `--fast` sets up Octopus for fast variant calling.

```shell
$ octopus -R ref.fa -I reads.bam --fast
```

**Warning** This command may reduce calling accuracy.

### `--very-fast`

Command `--very-fast` sets up Octopus for very fast variant calling.

```shell
$ octopus -R ref.fa -I reads.bam --very-fast

```

**Warning** This command may reduce calling accuracy

### `--data-profile`

Option `--data-profile` is used to generate a profile of the input read data, which may be used to make new [sequence error models](guides/errorModels.md).

```shell
$ octopus -R ref.fa -I reads.bam -o calls.bcf --data-profile reads.profile.csv

```

### `--bad-region-tolerance`

Option `--bad-region-tolerance` specifies the user tolerance for regions that may be 'uncallable' (e.g. due to mapping errors) and slow down calling. The possible arguments are:

* `LOW` Low tolerance to bad regions.
* `NORMAL` Default tolerance to bad regions.
* `HIGH` High tolerance to bad regions.
* `UNLIMITED` Turn off bad region detection.

```shell
$ octopus -R ref.fa -I reads.bam --bad-region-tolerance UNLIMITED

```

## Read pre-processing options

### `--disable-read-preprocessing`

Command `--disable-read-preprocessing` can be used to disable all optional read preprocessing - all viable raw input alignments will be used.

```shell
$ octopus -R ref.fa -I reads.bam --disable-read-preprocessing
```

### `--max-base-quality`

Option `--max-base-quality` caps the base quality of all input reads.

```shell
$ octopus -R ref.fa -I reads.bam --max-base-quality 40
```

### `--mask-low-quality-tails`

Option `--mask-low-quality-tails` is used to mask (assign base quality zero) read tail bases that have low base quality scores. The value provided to the option is the threshold used to define a 'low' quality score.

```shell
$ octopus -R ref.fa -I reads.bam --mask-low-quality-tails 5 // use 5 as minimum quality score
$ octopus -R ref.fa -I reads.bam --mask-low-quality-tails // use implicit 'low' score (see -h)
```

### `--mask-tails`

Option `--mask-tails` is used to unconditionally mask (assign base quality zero) read tail bases. The value provided to the option is the number of bases to mask.

```shell
$ octopus -R ref.fa -I reads.bam --mask-tails 3 // mask the last 3 bases of all reads
```

### `--mask-soft-clipped-bases`

Command `--mask-soft-clipped-bases` is used to mask (assign base quality zero) all read bases that are soft clipped in the input alignments.

```shell
$ octopus -R ref.fa -I reads.bam --mask-soft-clipped-bases
```

### `--soft-clip-mask-threshold`

Option `--soft-clip-mask-threshold` makes the option `--soft-clip-masking` conditional on the base quality score; only bases below the value given to the option will be masked.

```shell
$ octopus -R ref.fa -I reads.bam --soft-clip-masking yes --soft-clip-mask-threshold 10
```

### `--mask-soft-clipped-boundary-bases`

Option `--mask-soft-clipped-boundary-bases` makes the option `--soft-clip-masking` mask addition bases adjacent to soft clipped bases. The value provided to the option is the number of additional bases to mask.

```shell
$ octopus -R ref.fa -I reads.bam --soft-clip-masking yes --mask-soft-clipped-boundary-bases 5
```

### `--mask-inverted-soft-clipping`

Command `--mask-inverted-soft-clipping` is used to mask (assign base quality zero) all read bases that are soft clipped **and** are inverted copies of nearby sequence.

```shell
$ octopus -R ref.fa -I reads.bam --mask-inverted-soft-clipping
```

### `--mask-3prime-shifted-soft-clipped-heads`

Command `--mask-3prime-shifted-soft-clipped-heads` is used to mask (assign base quality zero) all read bases that are soft clipped **and** are shifted copies of nearby 3' sequence.

```shell
$ octopus -R ref.fa -I reads.bam --mask-3prime-shifted-soft-clipped-heads
```

### `--disable-adapter-masking`

Command `--disable-adapter-masking` is used to mask (assign base quality zero) read bases that are considered to be adapter contamination.

```shell
$ octopus -R ref.fa -I reads.bam --disable-adapter-masking
```

**Notes**

* The algorithm used to detect adapter contamination depends only on the input alignment mapping information; no library of adapter sequences are used.

### `--disable-overlap-masking`

Command `--disable-overlap-masking` is used to mask (assign base quality zero) read bases of read templates that contain overlapping segments, in order to remove non-independent base observations. 

```shell
$ octopus -R ref.fa -I reads.bam --disable-overlap-masking
```

**Notes**

* All but of the overlapping read bases are masked, leaving one base untouched. If two segments are overlapping, then half of the 5' bases of each segment are masked.

### `--split-long-reads`

Command `--split-long-reads` stipulates that reads longer than [`--max-read-length`](#option---max-read-length) should be split into smaller linked reads.

```shell
$ octopus -R ref.fa -I reads.bam --max-read-length 200 --split-long-reads
```

### `--consider-unmapped-reads`

Command `--consider-unmapped-reads` turns off the read filter that removes reads marked as unmapped in the input alignments.

```shell
$ octopus -R ref.fa -I reads.bam --consider-unmapped-reads
```

### `--min-mapping-quality`

Option `--min-mapping-quality` specifies the minimum mapping quality that reads must have to be considered; reads with mapping quality below this will be filtered and not used for analysis.

```shell
$ octopus -R ref.fa -I reads.bam --min-mapping-quality 10
```

### `--good-base-quality`

Option `--good-base-quality` defines the minimum quality of a 'good' base for the options `--min-good-base-fraction` and `--min-good-bases`.

```shell
$ octopus -R ref.fa -I reads.bam --good-base-quality 10
```

### `--min-good-base-fraction`

Option `--min-good-base-fraction` specifies the fraction of 'good' (see `--good-base-quality`) base qualities a read must have in order to be considered; reads with a fraction of 'good' base qualities less than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --min-good-base-fraction 0.3
```

### `--min-good-bases`

Option `--min-good-bases` specifies the number of 'good' (see `--good-base-quality`) base qualities a read must have in order to be considered; reads with 'good' base qualities less than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --min-good-bases 50
```

### `--allow-qc-fails`

Command `--allow-qc-fails` permits reads marked as QC fail in the input alignments to be considered, otherwise they will be filtered.

```shell
$ octopus -R ref.fa -I reads.bam --allow-qc-fails
```

### `--min-read-length`

Option `--min-read-length` specifies the minimum length (number of sequence bases) of reads to be considered; reads with length less than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --min-read-length 50
```

### `--max-read-length`

Option `--max-read-length` specifies the maximum length (number of sequence bases) of reads to be considered; reads with length greater than this will be filtered and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --max-read-length 500
```

### `--allow-marked-duplicates`

Command `--allow-marked-duplicates` disables the filter that removes reads marked as duplicate in the input alignments. The default behaviour is to remove reads marked duplicate.

```shell
$ octopus -R ref.fa -I reads.bam --allow-marked-duplicates
```

### `--allow-octopus-duplicates`

Command `--allow-octopus-duplicates` disables the filter that removes reads that Octopus considers to be duplicates. The default behaviour is to remove duplicate reads.

```shell
$ octopus -R ref.fa -I reads.bam --allow-octopus-duplicates
```

### `--duplicate-read-detection-policy`

Option `--duplicate-read-detection-policy` specifies approach to use for detecting [duplicate reads](guides/preprocessing.md). Possible arguments are:

* `RELAXED` Require 5' mapping co-ordinate matches **and** identical cigar strings.
* `AGGRESSIVE` Require 5' mapping co-ordinate matches only.

```shell
$ octopus -R ref.fa -I reads.bam --duplicate-read-detection-policy AGGRESSIVE
```

### `--allow-secondary-alignments`

Command `--allow-secondary-alignments` disables the filter that removes reads that are marked as secondary alignments in the input alignments. The default behaviour is to remove reads marked secondary.

```shell
$ octopus -R ref.fa -I reads.bam --allow-secondary-alignments
```

### `--allow-supplementary-alignments`

Command `--allow-supplementary-alignments` disables the filter that removes reads that are marked as supplementary alignments in the input alignments. The default behaviour is to remove reads marked supplementary.

```shell
$ octopus -R ref.fa -I reads.bam --allow-supplementary-alignments
```

### `--max-decoy-supplementary-alignment-mapping-quality`

Option `--max-decoy-supplementary-alignment-mapping-quality` removes any reads with supplementary alignments (i.e. `SA` SAM tag) to decoy contigs with mapping quality greater than the specified value.

```shell
$ octopus -R ref.fa -I reads.bam --max-decoy-supplementary-alignment-mapping-quality 20
```

### `--max-unplaced-supplementary-alignment-mapping-quality`

Command `--max-unplaced-supplementary-alignment-mapping-quality removes reads with supplementary alignments (i.e. `SA` SAM tag) to unplaced contigs with mapping quality greater than the specified value.

```shell
$ octopus -R ref.fa -I reads.bam --max-unplaced-supplementary-alignment-mapping-quality 20
```

### `--max-unlocalized-supplementary-alignment-mapping-quality`

Command `--max-unlocalized-supplementary-alignment-mapping-quality` removes reads with supplementary alignments (i.e. `SA` SAM tag) to unlocalized contigs with mapping quality greater than the specified value.

```shell
$ octopus -R ref.fa -I reads.bam --max-unlocalized-supplementary-alignment-mapping-quality 20
```

### `--no-reads-with-unmapped-segments`

Command `--no-reads-with-unmapped-segments` removes reads that are marked as having unmapped segments in the input alignments. The default behaviour is to allow reads with unmapped segments.

```shell
$ octopus -R ref.fa -I reads.bam --no-reads-with-unmapped-segments
```

### `--no-reads-with-distant-segments`

Command `--no-reads-with-distant-segments` removes reads that have segments mapped to different contigs. The default behaviour is to allow reads with distant segments.

```shell
$ octopus -R ref.fa -I reads.bam --no-reads-with-distant-segments
```

### `--no-adapter-contaminated-reads`

Command `--no-adapter-contaminated-reads` removes reads that are considered to have adapter contamination (i.e. sequence bases from the adapter). The default behaviour is to allow reads with adapter contamination.

```shell
$ octopus -R ref.fa -I reads.bam --no-adapter-contaminated-reads
```

### `--disable-downsampling`

Command `--disable-downsampling` turns off downsampling.

```shell
$ octopus -R ref.fa -I reads.bam --disable-downsampling
```

### `--downsample-above`

Option `--downsample-above` specifies the read depth required to mark a position as a candidate for the downsampler.

```shell
$ octopus -R ref.fa -I reads.bam --downsample-above 5000 # depths up to 5000x are permitted
```

### `--downsample-target`

Option `--downsample-target` specifies the target read depth for the downsampler for all candidate sites. Reads will be removed from the input alignments until all positions have read depth not greater than this.

```shell
$ octopus -R ref.fa -I reads.bam --downsample-target 100 # downsample to 100x
```

### `--use-same-read-profile-for-all-samples`

Command `--use-same-read-profile-for-all-samples` specifies that the same input read profile should be used for all samples, rather than generating one for each sample. This essentially means that the same read distribution is assumed for all samples. 

```shell
$ octopus -R ref.fa -I reads.bam --use-same-read-profile-for-all-samples
```

## Variant discovery options

### `--variant-discovery-mode`

Option `--variant-discovery-mode` specifies the thresholds used for candidate variant discovery, which affects the overall sensitivity of the generators. Possible values are:

* `ILLUMINA` The default mode, for Illumina quality reads.
* `PACBIO` For PacBio quality reads - requires more observations to propose a candidate, particular indel candidates.

```shell
$ octopus -R ref.fa -I reads.bam --variant-discovery-mode PACBIO
```

### `--disable-denovo-variant-discovery`

Command `--disable-denovo-variant-discovery` disables the pileup candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --disable-denovo-variant-discovery
```

### `--disable-repeat-candidate-generator`

Command `--disable-repeat-candidate-generator` disables the tandem repeat candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --disable-repeat-candidate-generator
```

### `--disable-assembly-candidate-generator`

Command `--disable-assembly-candidate-generator` disables the local *de novo* assembly candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --disable-assembly-candidate-generator
```

### `--source-candidates`

Option `--source-candidates` (short `-c`) accepts a list of VCF files, the contents of which will be added to the candidate variant set.

```shell
$ octopus -R ref.fa -I reads.bam --source-candidates calls.vcf.gz
```

**Notes**

* The option accepts files in the `.vcf`, `.vcf.gz`, and `.bcf` formats. Files in the `.vcf.gz` and `.bcf` format must be indexed. For larger

### `--source-candidates-file`

Option `--source-candidates-file` can be used to provide one or more files containing lists (one per line) of files in the `.vcf`, `.vcf.gz`, and `.bcf` formats (see notes on `--source-candidates`).

```shell
$ octopus -R ref.fa -I reads.bam --source-candidates-file vcfs.txt
```

where `vcfs.txt` may contain, for example

```text
/path/to/vcfs/callsA.vcf.gz
/path/to/vcfs/callsB.bcf
```

### `--min-source-candidate-quality`

Option `--min-source-candidate-quality` specifies the minimum `QUAL` score for a user provided source variant (see `--source-candidates` and `--source-candidates-file`) to be added to the final candidate variant list; any records with `QUAL` less than this will not be considered.

```shell
$ octopus -R ref.fa -I reads.bam -c calls.vcf.gz --min-source-candidate-quality 10
```

### `--use-filtered-source-candidates`

Command `--use-filtered-source-candidates` specifies allows variants in the user-provided candidate variants (see `--source-candidates` and `--source-candidates-file`) that are marked as filtered (according to the `FILTER` column) should be added to the final candidate variant list. The default behaviour is to remove filtered variants.

```shell
$ octopus -R ref.fa -I reads.bam -c calls.vcf.gz --use-filtered-source-candidates
```

### `--min-pileup-base-quality`

Option `--min-pileup-base-quality` specifies the minimum base quality that a SNV in the input alignments must have in ordered to be considered by the pileup candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --min-pileup-base-quality 10
```

### `--min-supporting-reads`

Option `--min-supporting-reads` specifies the minimum number of supporting reads a variant must have in the input alignments to be included in the candidate variant list from the pileup candidate variant generator.

```shell
$ octopus -R ref.fa -I reads.bam --min-supporting-reads 3
```

### `--allow-pileup-candidates-from-likely-misaligned-reads`

Command `--allow-pileup-candidates-from-likely-misaligned-reads` stops the pileup candidate generator filtering out candidates that are considered likely to originate from read misalignment, which can otherwise result in many false candidates.

```shell
$ octopus -R ref.fa -I reads.bam --allow-pileup-candidates-from-likely-misaligned-reads
```

### `--max-variant-size`

Option `--max-variant-size` specifies the maximum size (w.r.t reference region) a candidate variant can have to be considered by the calling algorithm; candidate variants with size greater than this will be removed and not considered.

```shell
$ octopus -R ref.fa -I reads.bam --max-variant-size 500
```

### `--kmer-sizes`

Option `--kmer-sizes` specifies the default kmer sizes to try for local *de novo* assembly. Assembly graphs will be constructed for each kmer size in all active assembler regions, and the union of candidate variants from each assembler graph used for the final candidate set.

```shell
$ octopus -R ref.fa -I reads.bam --kmer-sizes 5 10 15 20 25 30 35 40 45 50
```

### `--num-fallback-kmers`

Option `--num-fallback-kmers` specifies the number of additional kmer sizes to try in the case that none of the default (i.e. kmer sizes specified in `--kmer-sizes`) is able to construct a valid assembly graph (often the case for smaller kmer sizes).

```shell
$ octopus -R ref.fa -I reads.bam --num-fallback-kmers 20
```

### `--fallback-kmer-gap`

Option `--fallback-kmer-gap` speifies the increment between fallback kmer sizes (see `--num-fallback-kmers`)

```shell
$ octopus -R ref.fa -I reads.bam --fallback-kmer-gap 5
$ octopus -R ref.fa -I reads.bam --kmer-sizes 5 10 --num-fallback-kmers 2 --fallback-kmer-gap 5 # always try kmer sizes 5 and 10, and then try 15 and 20 if 5 and 10 don't work
```

### `--max-region-to-assemble`

Option `--max-region-to-assemble` specifies the maximum reference region size that can be used to construct a local *de novo* assembly graph. Larger values enable detection of larger variants (e.g. large deletions), but increase the likelihood that small kmer sizes will result in an invalid assembly graph, and therefore decrease sensitivity for smaller variation.

```shell
$ octopus -R ref.fa -I reads.bam --max-region-to-assemble 1000
```

### `--max-assemble-region-overlap`

Option `--max-assemble-region-overlap` specifies the maximum overlap between reference regions used to build local *de novo* assembly graphs. Larger overlaps result in more assembly graphs being constructed and may increase sensitivity for variation, at the expense of compute time. 

```shell
$ octopus -R ref.fa -I reads.bam --max-assemble-region-overlap 100
```

### `--assemble-all`

Command `--assemble-all` forces all reference regions to be used for local *de novo* assembly, rather than only regions considered likely to contain variation.

```shell
$ octopus -R ref.fa -I reads.bam --assemble-all
```

**Warning** This command may result in substantially longer runtimes. 

### `--assembler-mask-base-quality`

Option `--assembler-mask-base-quality` specifies the minimum base quality an aligned base should have to avoid being 'masked' (i.e. converted to reference) before being inserted into the local *de novo* assembly graph. Higher values reduce sensitivity to noise and increase the likelihood of finding bubbles in graphs with larger kmer sizes, as the cost of decreased sensitivity at lower kmer sizes. 

```shell
$ octopus -R ref.fa -I reads.bam --assembler-mask-base-quality 20
```

### `--allow-cycles`

Command `--allow-cycles` forces the assembler to consider assembly graphs containing non-reference cycles, which usually result in false candidates.

```shell
$ octopus -R ref.fa -I reads.bam --allow-cycles
```

### `--min-kmer-prune`

Option `--min-kmer-prune` specifies the minimum number of kmer observations that must be present in the local *de novo* assembly graph for the kmer to stay in the graph before variant extraction. Lower values increase sensitivity but lower specificity.

```shell
$ octopus -R ref.fa -I reads.bam --min-kmer-prune 3
```

### `--max-bubbles`

Option `--max-bubbles` specifies the maximum number of bubbles in the final local *de novo* assembly graph that may be explored for candidate variant generation. Higher values increase sensitivity but reduce specificity.

```shell
$ octopus -R ref.fa -I reads.bam --max-bubbles 100
```

### `--min-bubble-score`

Option `--min-bubble-score` specifies the minimum score a bubble explored in a local *de novo* assembly graph must have for the bubble to be extracted as a candidate variant. The 'bubble score' is the average of all kmer observations along the bubble edge, scaled by the probability the kmer observations along the bubble edge originate from a single read strand. The bubble score can be viewed as a proxy for the number of 'good' read observations for the bubble (i.e. variant), and therefore higher values increase specificity but reduce sensitivity.

```shell
$ octopus -R ref.fa -I reads.bam --min-bubble-score 5
```

### `--min-candidate-credible-vaf-probability`

Option `--min-candidate-credible-vaf-probability` sets the minimum probability mass above `--min-credible-somatic-frequency` required to 'discover' a variant when using the `cancer` calling model. Smaller values increase the number of candidate variants generated, potentially improving sensitivity but also increasing computational complexity.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-candidate-credible-vaf-probability 0.5
```

## Haplotype generation options

### `--max-haplotypes`

Option `--max-haplotypes` specifies the maximum number of haplotypes that are considered by the calling model. It also specifies the target number of haplotypes that the haplotype generator should produce on each iteration of the algorithm. If the haplotype generator is unable to satisfy the request (i.e. produces a greater number of haplotypes), then the set of haplotypes is reduced to this size by removing haplotypes considered unlikely using several likelihood based statistics.

Increasing `--max-haplotypes` reduces the chance that a true haplotype is incorrectly filtered before evaluation by the calling model. It also increases the number of candidate variants that may be considered on each iteration of the algorithm, potentially improving calling accuracy. However, increasing this value will usually increase runtimes - sometimes substantially. 

```shell
$ octopus -R ref.fa -I reads.bam --max-haplotypes 500
```

### `--haplotype-holdout-threshold`

Option `--haplotype-holdout-threshold` specifies the number of haplotypes that the haplotype generator can produce before some active candidate variants are added to the holdout stack. The value must not be less than `--max-haplotypes`.

```shell
$ octopus -R ref.fa -I reads.bam --haplotype-holdout-threshold 1000
```

### `--haplotype-overflow`

Option `--haplotype-overflow` specifies the number of haplotypes that the haplotype generator can produce before the current active region must be skipped. The value must not be less than `--haplotype-holdout-threshold`.

```shell
$ octopus -R ref.fa -I reads.bam --haplotype-overflow 5000
```

### `--max-holdout-depth`

Option `--max-holdout-depth` specifies the maximum size of the holdout stack.

```shell
$ octopus -R ref.fa -I reads.bam --max-holdout-depth 0 # no holdouts
```

### `--extension-level`

Option `--extension-level` specifies the condition for extending the active haplotype tree with novel alleles. The possible values are:

* `MINIMAL` Only include novel alleles overlapping with overlapping reads. 
* `NORMAL` Include novel alleles with reads overlapping the rightmost included allele.
* `AGGRESSIVE` No conditions on extension other than the number of alleles.

More aggressive extension levels result in larger haplotype blocks, which may improve phase lengths and accuracy. However, aggressive extension increases the possibility of including novel alleles that cannot be phased with active alleles which will increase compute time.

```shell
$ octopus -R ref.fa -I reads.bam --extension-level AGGRESSIVE
```

### `--lagging-level`

Option `--lagging-level` specifies the extent to which active alleles remain active in the next algorithm iteration, increasing the length of haplotype blocks.

* `NONE` Disable lagging; each iteration of the algorithm will evaluate a novel set of alleles.
* `NORMAL` Consider previous active alleles if there are reads overlapping them with the next active allele set.
* `AGGRESSIVE` Consider previous active alleles if there are overlapping reads that span them, and the next active alleles.

```shell
$ octopus -R ref.fa -I reads.bam --lagging-level AGGRESSIVE
```

### `--backtrack-level`

Option `--backtrack-level` specifies the extent to which 

* `NONE` Disable all backtracking.
* `NORMAL` Enables backtracking.
* `AGGRESSIVE` Currently the same as `NORMAL`.

```shell
$ octopus -R ref.fa -I reads.bam --backtrack-level NORMAL
```

### `--min-protected-haplotype-posterior`

Option `--min-protected-haplotype-posterior` specifies the minimum posterior probability that a haplotype is present in the samples (according to the calling model) for the haplotype to avoid being removed from consideration; haplotypes with posterior probability less than this may be filtered. Increasing the value of this option results in a greater number of haplotypes being filtered, allowing the haplotype tree to grow to include more candidate alleles, and reducing computational complexity. However, larger values also increase the chance that a true haplotype is incorrectly discarded. 

```shell
$ octopus -R ref.fa -I reads.bam --min-protected-haplotype-posterior 1e-5
```

### `--dont-protect-reference-haplotype`

Command `--dont-protect-reference-haplotype` disables protection of the reference haplotype during haplotype filtering, and therefore ensures that the reference haplotype is always considered by the calling model.

```shell
$ octopus -R ref.fa -I reads.bam --dont-protect-reference-haplotype
```

### `--bad-region-tolerance`

Option `--bad-region-tolerance` specifies the 'tolerance' for spending time calling variants in regions that are unlikely to be callable (e.g. due to low complexity sequence). Such regions tend to be computationally difficult and skipping them can save a lot of wasted computation time. However, identification of such regions is not perfect and could result in skipping regions with real variation. Possible argument values are:

* `LOW` Skip region that show any signs of being uncallable.
* `NORMAL` Skip region that show reasonable signs of being uncallable.
* `HIGH` Skip region that show strong signs of being uncallable.

```shell
$ octopus -R ref.fa -I reads.bam --bad-region-tolerance LOW
```

## Common variant calling options

### `--caller`

Option `--caller` (short '-C') specifies the calling model to be used. The option must only be set if the calling model is not automatically determined from other options.

```shell
$ octopus -R ref.fa -I reads.bam --caller cancer # e.g. for tumour-only
```

### `--organism-ploidy`

Option `--organism-ploidy` (short '-P') specifies the default ploidy of all input samples. All contigs will be assumed to have this ploidy unless specified otherwise in `--contig-ploidies`.

```shell
$ octopus -R ref.fa -I reads.bam --organism-ploidy 3 # triploid calling
```

### `--contig-ploidies`

Option `--contig-ploidies` (short `-p`) can be used to specify the ploidy of contigs (format `contig=ploidy`) that are not the same as the `--organism-ploidy`, or the ploidy of individual sample contigs (format `sample:contig=ploidy`).

```shell
$ octopus -R ref.fa -I reads.bam --contig-ploidies X=1 # Contig "X" has ploidy 1 for all samples
$ octopus -R ref.fa -I reads.bam --contig-ploidies HG002:X=1 # Contig "X" has ploidy 1 for sample "HG002"
```

### `--contig-ploidies-file`

Option `--contig-ploidies-file` can be used to provide a file specifying contig or sample contig ploidies (see `--contig-ploidies`), one per line.  

```shell
$ octopus -R ref.fa -I reads.bam --contig-ploidies-file ploidies.txt
```

where `ploidies.txt` could contain

```text
X=1
HG002:X=1
```

### `--min-variant-posterior`

Option `--min-variant-posterior` specifies the minimum posterior probability (Phred scale) required to call a variant. Candidate variants with posterior probability less than this will not be reported in the final call set.

```shell
$ octopus -R ref.fa -I reads.bam --min-variant-posterior 0.5
```

**Notes**

* For calling models with more than one class of variation, this option refers to the calling models default variant class. For example, in the `cancer` calling model there are both germline and somatic variants, and this option refers to germline variants (see the option `--min-somatic-posterior` for somatic variants).

### `--refcall`

Option `--refcall` can be used to enable reference calling, meaning Octopus will generate a 'gVCF' format output. There are two possible arguments for this option:

* `BLOCKED` Merge adjacent called reference positions with similar quality (see `--refcall-block-merge-quality`) a single gVCF record.
* `POSITIONAL` Emit a gVCF record for all positions.

If the option is specified but no argument is provide, `BLOCKED` is assumed.

```shell
$ octopus -R ref.fa -I reads.bam --refcall BLOCKED
$ octopus -R ref.fa -I reads.bam --refcall # same as above
$ octopus -R ref.fa -I reads.bam --refcall POSITIONAL
```

### `--refcall-block-merge-quality`

Option `--refcall-block-merge-quality` specifies the quality (Phred scale) threshold for merging adjacent called reference positions when `BLOCKED` refcalls are requested; adjacent reference positions with an absolute quality difference less or equal than this will be merged into a block.

```shell
$ octopus -R ref.fa -I reads.bam --refcall --refcall-block-merge-quality 20
```

### `--min-refcall-posterior`

Option `--min-refcall-posterior` specifies the minimum posterior probability required to call a position as homozygous reference; positions with posterior probability less than this will not be reported in the output gVCF.

```shell
$ octopus -R ref.fa -I reads.bam --min-refcall-posterior 10
```

### `--max-refcall-posterior`

Option `--max-refcall-posterior` caps the `QUAL` of all reference calls, which may result in larger reference blocks and smaller gVCF file sizes.

```shell
$ octopus -R ref.fa -I reads.bam --max-refcall-posterior 100
```

### `--snp-heterozygosity`

Option `--snp-heterozygosity` specifies the SNV heterozygosity parameter for the Coalesent mutation model, used to assign prior probabilities.

```shell
$ octopus -R ref.fa -I reads.bam --snp-heterozygosity 0.1
```

### `--snp-heterozygosity-stdev`

Option `--snp-heterozygosity-stdev` specifies the standard deviation of the SNV heterozygosity parameter (see `--snp-heterozygosity`).

```shell
$ octopus -R ref.fa -I reads.bam --snp-heterozygosity-stdev 0.1
```

### `--indel-heterozygosity`

Option `--indel-heterozygosity` specifies the INDEL heterozygosity parameter for the Coalesent mutation model, used to assign prior probabilities.

```shell
$ octopus -R ref.fa -I reads.bam --indel-heterozygosity 0.001
```

### `--max-genotypes`

Option `--max-genotypes` specifies the maximum number of candidate genotypes that must be evaluated by the calling model. If there are more possible candidate genotypes than this value then the algorithm may decide to remove some candidate genotypes using heuristics. The number of candidate genotypes to evaluate can have a substantial impact on runtime for some calling models (e.g. `cancer`).

```shell
$ octopus -R ref.fa -I reads.bam --max-genotypes 1000
```

### `--max-genotype-combinations`

Option `--max-genotype-combinations` specifies the maximum number of candidate joint genotype vectors that must be evaluated by the calling model. If there are more possible candidate joint genotype vectors than this value the algorithm may decide to remove some candidate genotype vectors using heuristics. The number of candidate joint genotype vectors to evaluate can have a substantial impact on runtime. 

```shell
$ octopus -R ref.fa -I reads.bam --max-genotype-combinations 100000
```

**Notes**

* This option is only used by calling models that consider joint genotypes (i.e. `population` and `trio`)

### `--use-uniform-genotype-priors`

Command `--use-uniform-genotype-priors` indicates that the uniform genotype prior model should be used for calling genotypes and variants.

```shell
$ octopus -R ref.fa -I reads.bam --use-uniform-genotype-priors
```

### `--use-independent-genotype-priors`

Command `--use-independent-genotype-priors` indicates that an independent genotype prior model should be used for evaluating joint genotype vectors.

```shell
$ octopus -R ref.fa -I reads.bam --use-independent-genotype-priors
```

### `--model-posterior`

Option `--model-posterior` enables or disables model posterior evaluation at each called variant site.

```shell
$ octopus -R ref.fa -I reads.bam --model-posterior yes
```

### `--disable-inactive-flank-scoring`

Command `--disable-inactive-flank-scoring` disables an additional step during haplotype likelihood calculation that attempts to correct low likelihoods caused by inactive variation in the flanking regions of the haplotype under evaluation.

```shell
$ octopus -R ref.fa -I reads.bam --disable-inactive-flank-scoring
```

### `--dont-model-mapping-quality`

Command `--dont-model-mapping-quality` disables consideration of mapping quality in the haplotype likelihood calculation. This can improve calling accuracy if read mapping qualities are well calibrated.

```shell
$ octopus -R ref.fa -I reads.bam --dont-model-mapping-quality
```

### `--sequence-error-model`

Option `--sequence-error-model` specifies the [sequence error model](guides/errorModels.md) to use.

```shell
$ octopus -R ref.fa -I reads.bam --sequence-error-model PCR.NOVASEQ # built-in error model
$ octopus -R ref.fa -I reads.bam --sequence-error-model /path/to/my/error.model # custom error model
```

**Notes**

* The same error model is used for all input reads. 

### `--max-indel-errors`

Option `--max-indel-errors` specifies the maximum number of indel errors in an individual read fragment that can be accuretely modelled by the haplotype likelihood model. Larger values usually require greater computational resources (determined by your systems available SIMD instructions).

```shell
$ octopus -R ref.fa -I reads.bam --max-indel-errors 32
```

### `--use-wide-hmm-scores`

Command `--use-wide-hmm-scores` sets the score variable computed by the pair HMM for haplotype likelihoods to 32 bits (the default is 16 bits). This can avoid score overflow in long noisy reads, but will slow does the computation.

```shell
$ octopus -R ref.fa -I reads.bam --use-wide-hmm-scores
```

### `--max-vb-seeds`

Option `--max-vb-seeds` specifies the maximum number of seeds that Variational Bayes models can use for posterior evaluation. Increasing the number of seeds increases the likelihood that a posterior mode will be identified, but results in more computation time.

```shell
$ octopus -R ref.fa -I reads.bam --max-vb-seeds 50
```

### `--read-linkage`

Option `--read-linkage` specifies how reads are linked in the input alignments. Read linkage information is used by the haplotype likelihood calculation and can improve calling accuracy and increase phase lengths.

* `NONE` Reads are not linked in any way.
* `PAIRED` Reads may be paired, with pairs having identical read names.
* `LINKED` Reads may be linked or paired, with linked reads having identical `BX` tags.

```shell
$ octopus -R ref.fa -I reads.bam --read-linkage LINKED
```

### `--min-phase-score`

Option `--min-phase-score` specifies the minimum phase score (`PQ` in VCF) required to emit adjacent variant calls in the same [phase set](guides/advanced/vcf.md#Haplotypes). Increasing this value results in less sites being phased, but reduces the phase false positive rate.

```shell
$ octopus -R ref.fa -I reads.bam --min-phase-score 5
```

### `--disable-early-phase-detection`

Command `--disable-early-phase-detection` prevents the phasing algorithm being applied to partially resolved haplotype blocks, which can lead to removal of complete phased segments from the head of the current haplotype block. This heuristic can prevent discontiguous phase blocks being resolved, which are more likely in some data (e.g. linked reads).

```shell
$ octopus -R ref.fa -I reads.bam --disable-early-phase-detection
```

## Cancer variant calling options

### `--normal-samples`

Option `--normal-samples` specifies which of the input samples are normal samples for tumour-normal paired analysis.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL
```

**Notes**

* Specifying this option will automatically activate the `cancer` calling model.

### `--max-somatic-haplotypes`

Option `--max-somatic-haplotypes` specifies the maximum number of unique haplotypes containing somatic variation that can be modelled by the somatic genotype model. If there are more true somatic haplotypes present in the input data than this value then the model will not accuretely fit the data and true somatic variants may not be called, however, larger values substantially increase the computational complexity of the model and potentially increase the false positive rate.

```shell
$ octopus -R ref.fa -I reads.bam --max-somatic-haplotypes 3
```

### `--somatic-snv-prior`

Option `--somatic-snv-prior` specifies the somatic SNV mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --somatic-snv-prior 1e-7
```

### `--somatic-indel-prior`

Option `--somatic-indel-prior` specifies the somatic INDEL mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --somatic-indel-prior 1e-7
```

### `--min-expected-somatic-frequency`

Option `--min-expected-somatic-frequency` specifies the minimum expected Variant Allele Frequency (VAF) for somatic mutations in the samples. This value is used by the `cancer` calling model as the lower bound for the VAF posterior marginalisation used to compute the posterior probability that a candidate variant is a somatic mutation. Decreasing this value increases sensitivity for somatic mutations with smaller VAFs, but also increases sensitivity to noise.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-expected-somatic-frequency 0.05
```

### `--min-credible-somatic-frequency`

Option `--min-credible-somatic-frequency` specifies the minimum credible Variant Allele Frequency (VAF) for somatic mutations in the samples. This value is used for candidate variant discovery, and also by the `cancer` calling model as the lower-bound on candidate somatic mutation VAF credible regions (i.e. if the credible region computed for a mutation contains VAFs less than this value then the mutation will not be called as somatic). Decreasing this option increases sensitivity for somatic mutations with lower VAFs, but also increases sensitivity to noise and increases computational complexity as more candidate variants will be generated.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-credible-somatic-frequency 0.001
```

**Notes**

* If the value provided to this option is greater than the value specified by `--min-expected-somatic-frequency`, then this value is used for both options.

### `--tumour-germline-concentration`

Option `--tumour-germline-concentration` sets the Dirichlet concentration parameter for the germline haplotypes of tumour samples. Larger values concentrate more prior probability mass on equal frequencies of germline haplotypes and also increase prior mass on low VAFs for somatic haplotypes. This can help the model correctly classify somatic variation when the normal sample is not informative (or not present), but also reduces sensitivity to somatic variation with larger VAFs.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --tumour-germline-concentration 5
```

### `--somatic-credible-mass`

Option `--somatic-credible-mass` specifies the probability mass to used to compute the credible interval for determining whether to call a variant as somatic (see also `--min-credible-somatic-frequency`). Larger values result in wider credible regions, increasing sensitivity for lower VAF somatic mutations, but also increase sensitivity to noise. 

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --somatic-credible-mass 0.99
```

### `--min-somatic-posterior`

Option `--min-somatic-posterior` specifies the minimum posterior probability (Phred scale) required to call a candidate variant as a somatic mutation.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --min-somatic-posterior 1
```

### `--normal-contamination-risk`

Option `--normal-contamination-risk` indicates the risk that the normal sample contains contamination from the tumour samples. There are two possible values:

* `LOW` The algorithm will not consider normal contamination when generating candidate genotypes.
* `HIGH` The algorithm will consider normal contamination when generating candidate genotypes.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --normal-contamination-risk HIGH
```

### `--somatics-only`

Command `--somatics-only` indicates that only variant sites called as somatic mutations (i.e. tagged `SOMATIC`) should appear in the final output.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --somatics-only
```

**Warning** Using this command will produce a VCF file that cannot be re-filtered using Octopus.

## Trio variant calling options

### `--maternal-sample`

Option `--maternal-sample` (short '-M') indicates which of the input samples is the mother of the proband. If this option is specified then `--paternal-sample` must also be specified.

```shell
$ octopus -R ref.fa -I reads.bam -M mother -F father
```

**Notes**

* Specifying this option will automatically activate the `trio` calling model.

### `--paternal-sample`

Option `--paternal-sample` (short '-F') indicates which of the input samples is the father of the proband. If this option is specified then `--maternal-sample` must also be specified.

```shell
$ octopus -R ref.fa -I reads.bam -M mother -F father
```

**Notes**

* Specifying this option will automatically activate the `trio` calling model.

### `--denovo-snv-prior`

Option `--denovo-snv-prior` specifies the *de novo* SNV mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovo-snv-prior 1e-7
```

### `--denovo-indel-prior`

Option `--denovo-indel-prior` specifies the *de novo* INDEL mutation prior probability for the samples under consideration.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovo-indel-prior 1e-8
```

### `--min-denovo-posterior`

Option `--min-denovo-posterior` specifies the minimum posterior probability (Phred scale) required to call a candidate variant as a *de novo* mutation.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --min-denovo-posterior 2
```

### `--denovos-only`

Command `--denovos-only` indicates that only variant sites called as *de novo* mutations (i.e. tagged `DENOVO`) should appear in the final output.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovos-only
```

**Warning** Using this command will produce a VCF file that cannot be re-filtered using Octopus.

## Polyclone variant calling options

### `--max-clones`

Option `--max-clones` specifies the maximum number of clones that can be modelled by the `polyclone` calling model. If there are more clones than this present in the input data then the model will not fit the data well and some true variation may not be called. However, larger values increase the computational complexity of the model and also increase sensitivity to noise.

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --max-clones 5
```

### `--min-clone-frequency`

Option `--min-clone-frequency` specifies the lower-bound Variant Allele Frequency (VAF) to use when computing the posterior probability for a variant. Smaller values increase sensitivity and the false positive rate.

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --min-clone-frequency 0.05
```

### `--clone-prior`

Option `--clone-prior` sets the prior probability for each new clone proposed by the model

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --clone-prior 0.1
```

### `--clone-concentration`

Option `--clone-concentration ` sets the concentration parameter for the symmetric Dirichlet distribution used to model clone frequencies

```shell
$ octopus -R ref.fa -I reads.bam --caller polyclone --clone-concentration 10
```

## Cell variant calling options

### `--max-copy-loss`

Option `--max-copy-loss` specifies the maximum number of germline haplotypes losses that can be considered by the model.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --max-copy-loss 1
```

### `--max-copy-gains`

Option `--max-copy-gains ` specifies the maximum number of haplotypes gains that can be considered by the model.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --max-copy-gains 1
```

### `--somatic-cnv-prior`

Option `--somatic-cnv-prior ` sets the prior probability of a loci having a copy-number change.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --somatic-cnv-prior 1e-10
```

### `--dropout-concentration`

Option `--dropout-concentration` sets the default Dirichlet concentration prior on haplotype frequencies. The higher the concentration parameter, the more probability mass around the centre of the distribution (0.5 for diploid). Lowering the concentration parameter means there’s more mass on extreme frequencies, which in turn means the model is less sensitive to dropout, but also less sensitive to real somatic variation. Setting to unity would put a uniform prior on frequencies.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --dropout-concentration 10
```

### `--sample-dropout-concentration`

Option `--sample-dropout-concentration` sets the Dirichlet concentration prior on haplotype frequencies for a specific sample, otherwise `--dropout-concentration` is used.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --sample-dropout-concentration NORMAL=100
```

### `--phylogeny-concentration`

Option `--phylogeny-concentration` sets the symmetric Dirichlet concentration prior on group mixture proportions in each phylogeny. A larger concentration parameter implies a more even distribution of samples across the tree.

```shell
$ octopus -R ref.fa -I reads.bam --caller cell --phylogeny-concentration 10
```

## Variant filtering options

### `--disable-call-filtering`

Command `--disable-call-filtering` disables variant call filtering.

```shell
$ octopus -R ref.fa -I reads.bam --disable-call-filtering
```

### `--filter-expression`

Option `--filter-expression` sets the threshold filter expression to use for filtering variants not tagged with `SOMATIC` or `DENOVO`.

```shell
$ octopus -R ref.fa -I reads.bam --filter-expression "QUAL < 10" # PASS calls with QUAL >= 10
```

### `--somatic-filter-expression`

Option `--somatic-filter-expression` sets the threshold filter expression to use for filtering `SOMATIC` variants.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL --somatic-filter-expression "QUAL < 20 | PP < 20"
```

### `--denovo-filter-expression`

Option `--denovo-filter-expression` sets the threshold filter expression to use for filtering `DENOVO` variants.

```shell
$ octopus -R ref.fa -I reads.bam --ped trio.ped --denovo-filter-expression "QUAL < 20 | PP < 20"
```

### `--refcall-filter-expression`

Option `--refcall-filter-expression` sets the threshold filter expression to use for filtering sites called homozygous reference.

```shell
$ octopus -R ref.fa -I reads.bam --refcall --refcall-filter-expression "QUAL < 10"
```

### `--use-preprocessed-reads-for-filtering`

Command `--use-preprocessed-reads-for-filtering` forces use of the same read pre-processing steps used for calling variants for filtering variants; otherwise all well-formed reads are used for filtering.

```shell
$ octopus -R ref.fa -I reads.bam --use-preprocessed-reads-for-filtering
```

### `--keep-unfiltered-calls`

Command `--keep-unfiltered-calls` requests that Octopus keep a copy of the VCF file produced before filtering is applied. The copy has `unfiltered` appended to the final output name.

```shell
$ octopus -R ref.fa -I reads.bam -o calls.vcf --keep-unfiltered-calls # writes calls.unfiltered.vcf
```

### `--annotations`

Option `--annotations` requests that the values of a sub-set of the measures used for filtering are reported in the final VCF output.

```shell
$ octopus -R ref.fa -I reads.bam --annotations SB # adds SB FORMAT field to each record
```

### `--filter-vcf`

Option `--filter-vcf` specifies an Octopus VCF file to filter. No calling is performed.

```shell
$ octopus -R ref.fa -I reads.bam --filter-vcf calls.bcf
```

### `--forest-model`

Option `--forest-model` enables [random forest variant filtering](guides/filtering/forest.md). The argument to the option is a ranger forest file.

```shell
$ octopus -R ref.fa -I reads.bam \
	--forest-model octopus/resources/forests/germline.v0.6.3-beta.forest
```

### `--somatic-forest-model`

Option `--somatic-forest-model` enables [random forest variant filtering](guides/filtering/forest.md) for somatic variants. The argument to the option is a ranger forest file.

```shell
$ octopus -R ref.fa -I reads.bam -N NORMAL \
	--forest-model octopus/resources/forests/germline.v0.6.3-beta.forest
	--somatic-forest-model octopus/resources/forests/somatic.v0.6.3-beta.forest
```

### `--min-forest-quality`

Option `--min-forest-quality` specifies the random forest minimum quality score (phred scale) required to PASS a variant call (`RFGQ_ALL`), and each samples genotype calls (`RFGQ`).

```shell
$ octopus -R ref.fa -I reads.bam --min-forest-quality 7
```

### `--use-germline-forest-for-somatic-normals`

Command `--use-germline-forest-for-somatic-normals` specifies that the forest model given to `--forest-model` should be used to score normal sample genotypes for in somatic records, rather than the forest model given to `--somatic-forest-model`.

```shell
$ octopus -R ref.fa -I normal.bam tumour.bam -N NORMAL --forest germline.forest --somatic-forest somatic.forest -o calls.vcf --use-germline-forest-for-somatic-normals
```