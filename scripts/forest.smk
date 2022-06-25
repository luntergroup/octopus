from pathlib import Path
import re
import pysam as ps
import csv
import subprocess as sp
import numpy as np
import random

default_prefix = "octopus"
default_germline_measures = "AC AD ADP AF AFB ARF BMQ BQ CC CRF DAD DAF DC DENOVO DP DPC ER ERS FRF GC GQ GQD ITV MC MF MHL MP MRC MQ MQ0 MQD PLN PP PPD QD QUAL REB RSB RTB SB SD SF STRL STRP VL".split()
default_somatic_measures = "AC AD ADP AF ARF BMQ BQ CC CRF DAD DAF DP DPC ER ERS FRF GC GQ GQD ITV NC MC MF MHL MP MRC MQ MQ0 MQD PLN PP PPD QD QUAL REB RSB RTB SB SD SF SHC SMQ SOMATIC STRL STRP VL".split()

known_truth_set_urls = {
    "GIAB": {
        "GRCh37": {
            "HG001": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
            },
            "HG002": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed"
            },
            "HG003": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG003_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG003_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG003_GRCh37_1_22_v4.1_draft_benchmark.bed"
            },
            "HG004": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG004_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG004_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG004_GRCh37_1_22_v4.1_draft_benchmark.bed"
            },
            "HG005": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf_noMetaSV.bed"
            }
        },
        "GRCh38": {
            "HG001": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
            },
            "HG002": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG002_GRCh38_1_22_v4.2_benchmark.bed"
            },
            "HG003": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG003_GRCh38_1_22_v4.2_benchmark.bed"
            },
            "HG004": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/data/AshkenazimTrio/analysis/NIST_v4.2_SmallVariantDraftBenchmark_07092020/HG004_GRCh38_1_22_v4.2_benchmark.bed"
            },
            "HG005": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.bed"
            }
        }
    }
}

for reference in ["GRCh37", "GRCh38"]:
    sample_synonyms = {"HG001": "NA12878", "HG002": "NA24385", "HG005": "NA24631"}
    for sample, synonym in sample_synonyms.items():
        known_truth_set_urls["GIAB"][reference][synonym] = known_truth_set_urls["GIAB"][reference][sample]

class TrainingExample:
    def __init__(self, data):
        self.kind = data["kind"] if "kind" in data else "germline"
        self.reference = Path(data["reference"]) if "reference" in data else None
        self.reads = data["reads"] if "reads" in data else None
        if type(self.reads) is not list:
            self.reads = [self.reads]
        self.reads = [Path(r) for r in self.reads]
        self.name = data["name"] if "name" in data else "+".join(bam.stem for bam in self.reads)
        self.regions = Path(data["regions"]) if "regions" in data else None
        self.truth = data["truth"] if "truth" in data else None
        self.confident = Path(data["confident_regions"]) if "confident_regions" in data else None
        self.tp_fraction = data["tp_fraction"] if "tp_fraction" in data else 1
        self.fp_fraction = data["fp_fraction"] if "fp_fraction" in data else 1
        self.threads = int(data["threads"]) if "threads" in data else 20
        self.options = data["options"] if "options" in data else []
        if type(self.options) is str:
            self.options = self.options.split(' ')

class ForestHyperparameters:
    def __init__(self, d):
        self.trees = d["trees"] if "trees" in d else 100
        self.min_samples_split = d["min_samples_split"] if "min_samples_split" in d else 1
        self.max_depth = d["max_depth"] if "max_depth" in d else 0

class TrainingOptions:
    def __init__(self, d):
        self.training_fraction = d["training_fraction"] if d is not None and "training_fraction" in d else None
        self.hyperparameters = ForestHyperparameters(d["hyperparameters"] if d is not None and "hyperparameters" in d else {})

prefix = default_prefix
measures = default_germline_measures

def load_config():
    outdir = config["dir"] if "dir" in config else ""
    global prefix
    if "prefix" in config:
        prefix = config["prefix"]

    examples = config['examples'] if "examples" in config else config
    if type(examples) is not list:
        examples = [examples]
    examples = [TrainingExample(example) for example in examples]
    
    global measures
    if "measures" in config:
        measures = config["measures"]
        if type(measures) is str:
            measures = measures.split(" ")
    elif all(example.kind == "somatic" for example in examples):
        measures = default_somatic_measures
    elif any(example.kind == "somatic" for example in examples):
        measures = list(set(default_germline_measures + default_somatic_measures))
    else:
        measures = default_germline_measures

    sdf_references, known_truths, given_truths = {}, {}, config["truths"] if "truths" in config else {}

    for label, details in given_truths.items():
        for detail, path in details.items():
            given_truths[label][detail] = Path(path)
        if not details["vcf"].exists():
            raise ValueError(str(details["vcf"]) + " does not exist")
        if not details["bed"].exists():
            raise ValueError(str(details["bed"]) + " does not exist")

    for example in examples:
        if type(example.truth) is OrderedDict:
            if example.confident is None:
                example.confident = {}
            for sample, truth_name in example.truth.items():
                if truth_name in given_truths:
                    example.truth[sample], example.confident[sample] = given_truths[truth_name]['vcf'], given_truths[truth_name]['bed']
                elif truth_name in known_truths:
                    example.truth[sample], example.confident[sample] = known_truths[example.truth]
                else:
                    library, reference_version, truth_sample = truth_name.split('//')
                    truth_prefix = Path("truth") / (truth_sample + "." + reference_version + "." + library + ".vcf")
                    vcf, bed = truth_prefix.with_suffix(".vcf.gz"), truth_prefix.with_suffix(".bed")
                    known_truths[truth_name] = vcf, bed
                    example.truth[sample], example.confident[sample] = vcf, bed
        else:
            if not Path(example.truth).exists():
                if example.truth in given_truths:
                    example.truth, example.confident = given_truths[example.truth]['vcf'], given_truths[example.truth]['bed']
                elif example.truth in known_truths:
                    example.truth, example.confident = known_truths[example.truth]
                else:
                    library, reference_version, truth_sample = example.truth.split('//')
                    truth_prefix = Path("truth") / (truth_sample + "." + reference_version + "." + library + ".vcf")
                    vcf, bed = truth_prefix.with_suffix(".vcf.gz"), truth_prefix.with_suffix(".bed")
                    known_truths[example.truth] = vcf, bed
                    example.truth, example.confident = vcf, bed
            else:
                if not example.confident.exists():
                    raise ValueError(example.confident + " does not exist")
            read_samples = [ps.AlignmentFile(example.reads[n]).header["RG"][0]["SM"] for n in range(len(example.reads))]
            example.truth, example.confident = {sample: example.truth for sample in read_samples}, {sample: example.confident for sample in read_samples}

    return {example.name: example for example in examples}, TrainingOptions(config['training'] if "training" in config else None)

examples, options = load_config()

rule all:
    input:
        prefix + ".forest"
localrules: all

rule samtools_faidx:
    input:
        "{name}.fa"
    output:
        "{name}.fa.fai"
    shell:
        "samtools faidx {input}"
localrules: samtools_faidx

def get_threads(example):
    if "--threads" in example.options:
        return int(example.options[example.options.index("--threads") + 1])
    else:
        return 1

rule octopus_call:
    input:
        reference = lambda wildcards: str(examples[wildcards.label].reference),
        reference_index = lambda wildcards: str(examples[wildcards.label].reference) + ".fai",
        bams = lambda wildcards: [str(bam) for bam in examples[wildcards.label].reads],
        bam_indices = lambda wildcards: [str(bam) + ".bai" for bam in examples[wildcards.label].reads],
        regions = lambda wildcards: str(examples[wildcards.label].regions)
    output:
        vcf="calls/{label}.octopus.raw.vcf.gz",
        vcf_index="calls/{label}.octopus.raw.vcf.gz.tbi"
    params:
        temp_dir = lambda wildcards: "calls/" + wildcards.label + "octopus.raw.temp",
        other = lambda wildcards: " ".join(examples[wildcards.label].options)
    threads: lambda wildcards: examples[wildcards.label].threads
    log:
        "logs/{label}.octopus.raw.log"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bams} \
         -t {input.regions} \
         -o {output.vcf} \
         --disable-call-filter \
         --threads {threads} \
         --temp-directory-prefix {params.temp_dir} \
         {params.other} \
        )2> {log}"

rule octopus_annotate:
    input:
        reference = lambda wildcards: str(examples[wildcards.label].reference),
        reference_index = lambda wildcards: str(examples[wildcards.label].reference) + ".fai",
        bams = lambda wildcards: [str(bam) for bam in examples[wildcards.label].reads],
        bam_indices = lambda wildcards: [str(bam) + ".bai" for bam in examples[wildcards.label].reads],
        regions = lambda wildcards: str(examples[wildcards.label].regions),
        vcf = rules.octopus_call.output.vcf,
        vcf_index = rules.octopus_call.output.vcf_index
    output:
        vcf = "calls/{label}.octopus.ann.vcf.gz",
        vcf_index = "calls/{label}.octopus.ann.vcf.gz.tbi"
    params:
        temp_dir = lambda wildcards: "calls/" + wildcards.label + "octopus.ann.temp",
        other = lambda wildcards: " ".join(examples[wildcards.label].options)
    threads: lambda wildcards: examples[wildcards.label].threads
    log:
        "logs/{label}.octopus.ann.log"
    shell:
        "(octopus \
         -R {input.reference} \
         -I {input.bams} \
         -t {input.regions} \
         -o {output.vcf} \
         --filter-vcf {input.vcf} \
         --disable-call-filter \
         --annotations all \
         --aggregate-annotations \
         --threads {threads} \
         --temp-directory-prefix {params.temp_dir} \
         {params.other} \
         )2> {log}"

rule rename_germline:
    input:
        vcf = rules.octopus_annotate.output.vcf,
        vcf_index = rules.octopus_annotate.output.vcf_index
    output:
        vcf = "calls/{label}.octopus.ann.germline.vcf.gz",
        vcf_index = "calls/{label}.octopus.ann.germline.vcf.gz.tbi"
    shell:
        "mv {input.vcf} {output.vcf} && mv {input.vcf_index} {output.vcf_index}"
localrules: rename_germline

rule filter_somatics:
    input:
        vcf = rules.octopus_annotate.output.vcf,
        vcf_index = rules.octopus_annotate.output.vcf_index
    output:
        vcf = "calls/{label}.octopus.ann.somatic.vcf.gz",
        vcf_index = "calls/{label}.octopus.ann.somatic.vcf.gz.tbi"
    shell:
        "bcftools view -i FORMAT/SOMATIC[*]=1 -Oz -o {output.vcf} {input.vcf} && tabix {output.vcf}"
localrules: filter_somatics

rule rtg_format:
    input:
        fa="{name}.fa",
        fai="{name}.fa.fai"
    output:
        directory("{name}.sdf")
    shell:
        "rtg format -o {output} {input.fa}"  
localrules: rtg_format

rule download_giab:
    output:
        vcf = "truth/{sample}.{reference}.GIAB.vcf.gz",
        vcf_index = "truth/{sample}.{reference}.GIAB.vcf.gz.tbi",
        bed = "truth/{sample}.{reference}.GIAB.bed"
    params:
        vcf = lambda wildcards: known_truth_set_urls["GIAB"][wildcards.reference][wildcards.sample]["vcf"],
        bed = lambda wildcards: known_truth_set_urls["GIAB"][wildcards.reference][wildcards.sample]["bed"]
    shell:
        """
        curl -o {output.vcf} {params.vcf}
        curl -o {output.vcf_index} {params.vcf}.tbi
        curl -o {output.bed} {params.bed}
        """
localrules: download_giab

def get_vcfeval_ploidy_option(example, sample):
    if "--organism-ploidy" in example.options:
        ploidy_idx = example.options.index("--organism-ploidy")
        if ploidy_idx >= 0:
            ploidy = example.options[ploidy_idx + 1]
            return "--sample-ploidy=" + str(ploidy)
    return ""

rule vcfeval:
    input:
        baseline_vcf = lambda wildcards: str(examples[wildcards.label].truth[wildcards.sample]),
        baseline_vcf_idx = lambda wildcards: str(examples[wildcards.label].truth[wildcards.sample]) + ".tbi",
        calls_vcf = "calls/{label}.octopus.ann.{kind}.vcf.gz",
        calls_vcf_idx = "calls/{label}.octopus.ann.{kind}.vcf.gz.tbi",
        reference = lambda wildcards: str(examples[wildcards.label].reference.with_suffix(".sdf"))
    output:
        directory("eval/{label}.octopus.ann.{sample}.{kind}.vcfeval")
    params:
        evaluation_regions = lambda wildcards: str(examples[wildcards.label].confident[wildcards.sample]),
        bed_regions = lambda wildcards: str(examples[wildcards.label].regions),
        squash_ploidy = lambda wildcards: "--squash-ploidy" if wildcards.kind == "somatic" else "",
        ploidy = lambda wildcards: get_vcfeval_ploidy_option(examples[wildcards.label], wildcards.sample),
        calls_sample = lambda wildcards: "ALT" if wildcards.kind == "somatic" else wildcards.sample,
        mode = "annotate"
    threads: 20
    shell:
        "bcftools view -h {input.baseline_vcf} | \
         tail -1 | awk '{{if($NF==\"INFO\"){{print \"ALT\"}}else{{print $NF}}}}' | \
         xargs -I{{}} \
         rtg vcfeval \
         -t {input.reference} \
         -b {input.baseline_vcf} \
         --bed-regions {params.bed_regions} \
         --evaluation-regions {params.evaluation_regions} \
         -c {input.calls_vcf} \
         --sample {{}},{params.calls_sample} \
         --all-records \
         --ref-overlap \
         --no-roc \
         --threads {threads} \
         -m {params.mode} \
         -o {output} \
         {params.squash_ploidy} \
         {params.ploidy} \
         "

def is_homref(vcf_rec, sample):
    return all(allele == vcf_rec.ref for allele in vcf_rec.samples[sample].alleles)

def has_homref_calls(vcf_filename, sample=None):
    vcf = ps.VariantFile(str(vcf_filename))
    if sample is None:
        assert len(vcf.header.samples) == 1
        sample = vcf.header.samples[0]
    for rec in vcf:
        if is_homref(rec, sample): return True
    return False

def complement_vcf(src_vcf_filename, tagret_vcf_filenames, dst_vcf_filename, regions_bed=None):
    cmd = ['bcftools', 'isec', '-C', str(src_vcf_filename)] + [str(f) for f in tagret_vcf_filenames] + ['-w1', '-Oz', '-o', str(dst_vcf_filename)]
    if regions_bed is not None:
        cmd += ["-R", regions_bed]
    sp.call(cmd)
    index_vcf(dst_vcf_filename)

def intersect_vcfs(src_vcf_filenames, dst_vcf_filename, regions_bed=None):
    cmd = ['bcftools', 'isec'] + [str(f) for f in src_vcf_filenames] + ['-n', str(len(src_vcf_filenames)), '-w1', '-Oz', '-o', str(dst_vcf_filename)]
    if regions_bed is not None:
        cmd += ["-R", regions_bed]
    sp.call(cmd)
    index_vcf(dst_vcf_filename)

def concat_vcfs(vcfs, out, remove_duplicates=True):
    assert len(vcfs) > 1
    cmd = ['bcftools', 'concat', '-a', '-Oz', '-o', str(out)]
    if remove_duplicates:
        cmd.append('-D')
    cmd += [str(vcf) for vcf in vcfs]
    sp.call(cmd)
    index_vcf(out)

def index_vcf(vcf_filename):
    sp.call(['tabix', '-f', str(vcf_filename)])

def remove_vcf_index(vcf_filename):
    vcf_index_filename = vcf_filename.with_suffix(vcf_filename.suffix + '.tbi')
    if vcf_index_filename.exists(): vcf_index_filename.unlink()

def remove_vcf(vcf_filename, remove_index=True):
    vcf_filename.unlink()
    if remove_index: remove_vcf_index(vcf_filename)

def rename_vcf(old_name, new_name, index=True):
    old_name.rename(new_name)
    if index:
        old_index_name = old_name.with_suffix(old_name.suffix + ".tbi")
        if old_index_name.exists():
            new_index_name = new_name.with_suffix(new_name.suffix + ".tbi")
            old_index_name.rename(new_index_name)

def classify_vcfeval_ignored_calls(vcfeval_dir, octopus_vcf, sample, regions=None):
    calls_vcf = ps.VariantFile(vcfeval_dir / "calls.vcf.gz")
    baseline_vcf = ps.VariantFile(vcfeval_dir / "baseline.vcf.gz")
    new_calls_vcf = ps.VariantFile(vcfeval_dir / "calls.homref.vcf.gz", 'wz', header=calls_vcf.header)
    for call in calls_vcf:
        if call.info["CALL"] == "IGN" and call.contig in baseline_vcf.header.contigs:
            if any(baseline.info["BASE"] == "FN" for baseline in baseline_vcf.fetch(call.contig, call.start, call.stop)):
                call.info["CALL"] = "FP"
            else:
                call.info["CALL"] = "TP"
        new_calls_vcf.write(call)
    calls_vcf.close()
    baseline_vcf.close()
    new_calls_vcf.close()
    rename_vcf(vcfeval_dir / "calls.homref.vcf.gz", vcfeval_dir / "calls.vcf.gz")

rule classify_homref_calls:
    input:
        eval_dir = "eval/{label}.octopus.ann.{sample}.{kind}.vcfeval",
        calls_vcf = "calls/{label}.octopus.ann.{kind}.vcf.gz"
    output:
        temporary("eval/{label}.octopus.ann.{sample}.{kind}.homref.done")
    params:
        regions = lambda wildcards: examples[wildcards.label].confident[wildcards.sample]
    run:
        classify_vcfeval_ignored_calls(Path(input[0]), Path(input[1]), wildcards.sample, params.regions)
        open(output[0], 'w')
localrules: classify_homref_calls

def get_annotation(field, rec, sample=None):
    if field == 'QUAL':
        return rec.qual
    elif field in rec.format:
        if sample is None:
            assert len(rec.samples) == 1
            res = rec.samples[list(rec.samples)[0]][field]
        else:
            res = rec.samples[sample][field]
        if type(res) == tuple:
            res = np.mean(res)
        return res
    elif field in rec.info:
        res = rec.info[field]
        if type(res) == tuple:
            res = np.mean(res)
        return res
    else:
        # Field must be a flag and not present
        return 0

def is_missing(x):
    return x is None or x == '.' or np.isnan(float(x))

def to_str(x):
    return str(x) if type(x) != bool else str(int(x))

def annotation_to_string(x, missing_value):
    return to_str(missing_value) if is_missing(x) else to_str(x)

def encode_classification(classifcation):
    if classifcation == "TP":
        return 1
    elif classifcation == "FP":
        return 0
    else:
        return 2

def make_ranger_data(vcf_filename, out_filename, sample, classifcation, measures, missing_value=-1, fraction=1):
    vcf = ps.VariantFile(vcf_filename)
    with out_filename.open(mode='w') as ranger_data:
        datawriter = csv.writer(ranger_data, delimiter=' ')
        encoded_classification = encode_classification(classifcation)
        for rec in vcf:
            if rec.info["CALL"] == classifcation and (fraction >= 1 or random.random() <= fraction):
                row = [annotation_to_string(get_annotation(measure, rec, sample=sample), missing_value) for measure in measures]
                row.append(encoded_classification)
                datawriter.writerow(row)

rule extract_annotations:
    input:
        rules.vcfeval.output,
        rules.classify_homref_calls.output
    output:
        "dat/{label}.octopus.ann.{sample}.{kind}.{classification}.dat"
    params:
        vcf = lambda _, input: Path(input[0]) / "calls.vcf.gz",
        measures = measures,
        fraction = lambda wildcards: examples[wildcards.label].tp_fraction if wildcards.kind.startswith("tp") else  examples[wildcards.label].fp_fraction
    wildcard_constraints:
        kind = "|".join([re.escape(s) for s in ["germline", "somatic"]]),
        classification = "|".join([re.escape(s) for s in ["TP", "FP_CA", "FP"]])
    run:
        make_ranger_data(params.vcf, Path(output[0]), wildcards.sample, wildcards.classification, params.measures, fraction=params.fraction) 
localrules: extract_annotations

rule cat_annotations:
    input:
        expand("dat/{{label}}.octopus.ann.{{sample}}.{{kind}}.{classification}.dat",
               classification=["TP", "FP_CA", "FP"])
    output:
        "dat/{label}.octopus.ann.{sample}.{kind}.dat"
    wildcard_constraints:
        kind = "|".join([re.escape(s) for s in ["germline", "somatic"]])
    shell:
        "cat {input} > {output}"
localrules: cat_annotations

def get_sample_kind(example, sample):
    if example.kind == "somatic" and "--normal-samples" in example.options:
        normal_sample_idx = example.options.index("--normal-samples")
        if example.options[normal_sample_idx + 1] == sample:
            return "germline"
    return example.kind

rule cat_data:
    input:
        ["dat/" + example.name + ".octopus.ann." + sample + "." + get_sample_kind(example, sample) + ".dat" \
         for _, example in examples.items() \
         for sample, _ in example.truth.items()]
    output:
        "dat/" + prefix + ".full.dat"
    shell:
        "cat {input} > {output}"
localrules: cat_data

def make_ranger_master_data_file(in_filename, out_filename, header, sample=None):
    with in_filename.open() as in_file, out_filename.open(mode='w') as out_file:
        out_file.write(header + '\n')
        for line in in_file:
            if sample is None or random.random() <= sample:
                out_file.write(line)

rule make_ranger_data:
    input:
        rules.cat_data.output
    output:
        "dat/" + prefix + ".training.dat"
    params:
        header = " ".join(measures + ['TP']),
        training_fraction = options.training_fraction
    run:
        make_ranger_master_data_file(Path(input[0]), Path(output[0]), params.header, sample=params.training_fraction)
localrules: make_ranger_data

rule install_ranger:
    output:
        "ranger/cpp_version/build/ranger"
    shell:
        """
        rm -rf ranger
        git clone https://github.com/dancooke/ranger
        mkdir ranger/cpp_version/build && cd "$_"
        cmake ..
        make
        """

localrules: install_ranger

rule ranger_train:
    input:
        ranger = rules.install_ranger.output,
        data = rules.make_ranger_data.output
    output:
        prefix + ".forest"
    params:
        trees = options.hyperparameters.trees,
        target_partition_size = options.hyperparameters.min_samples_split,
        max_depth = options.hyperparameters.max_depth,
        prefix = prefix
    threads: 100
    log:
        "logs/" + prefix + ".ranger.log"
    shell:
        "({input.ranger} \
         --file {input.data} \
         --depvarname TP \
         --probability \
         --nthreads {threads} \
         --outprefix {params.prefix} \
         --write \
         --impmeasure 1 \
         --ntree {params.trees} \
         --targetpartitionsize {params.target_partition_size} \
         --maxdepth {params.max_depth} \
         --verbose \
        )2> {log}"
