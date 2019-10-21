#!/usr/bin/env python3

import argparse
from os import makedirs, remove, pardir
from os.path import join, basename, exists, dirname, abspath
import subprocess as sp
import csv
import pysam as ps
import random
import numpy as np
import shutil
from sklearn.metrics import log_loss
import json
import urllib.request

script_dir = dirname(abspath(__file__))
default_octopus_bin = join(abspath(join(script_dir, pardir)), 'bin/octopus')

default_measures = "AC AD ADP AF ARF BQ CC CRF DAD DAF DP DPC ER ERS FRF GC GQ GQD NC MC MF MP MRC MQ MQ0 MQD PP PPD QD QUAL REFCALL REB RSB RTB SB SD SF SHC SMQ SOMATIC STR_LENGTH STR_PERIOD VL".split()

known_truth_set_urls = {
    "GIAB": {
        "GRCh37": {
            "NA12878": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
            },
            "HG001": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed"
            },
            "NA24385": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
            },
            "HG002": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
            },
            "NA24631": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf_noMetaSV.bed"
            },
            "HG005": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh37/HG005_GRCh37_highconf_CG-IllFB-IllGATKHC-Ion-SOLID_CHROM1-22_v.3.3.2_highconf_noMetaSV.bed"
            }
        },
        "GRCh38": {
            "NA12878": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
            },
            "HG001": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
            },
            "NA24385": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
            },
            "HG002": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed"
            },
            "NA24631": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.bed"
            },
            "HG005": {
                'vcf': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz",
                'vcf_idx': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.vcf.gz.tbi",
                'bed': "ftp://ftp-trace.ncbi.nlm.nih.gov//giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv3.3.2/GRCh38/HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf.bed"
            }
        }
    }
}

def get_octopus_version(octopus_bin):
    version_line = sp.check_output([octopus_bin, '--version']).decode("utf-8").split('\n')[0].split()
    return version_line[2], version_line[-1][:-1] if len(version_line) > 3 else None

def check_exists(paths):
    for path in paths: assert exists(path)

def check_exists_or_none(paths):
    for path in paths: assert path is None or exists(path)

class TrainingData:
    def __init__(self, data):
        self.reference = data["reference"] if "reference" in data else None
        self.sdf = data["reference_sdf"] if "reference_sdf" in data else None
        self.reads = data["reads"] if "reads" in data else None
        if type(self.reads) is not list:
            self.reads = [self.reads]
        self.octopus_vcf = data["octopus_vcf"] if "octopus_vcf" in data else None
        self.regions = data["calling_regions"] if "calling_regions" in data else None
        self.truth = data["truth"] if "truth" in data else None
        self.confident = data["confident_regions"] if "confident_regions" in data else None
        self.tp = data["tp_fraction"] if "tp_fraction" in data else 1
        self.fp = data["fp_fraction"] if "fp_fraction" in data else 1
        self.config = data["octopus_config"] if "octopus_config" in data else None

class TrainingOptions:
    def __init__(self, d):
        self.cross_validation_fraction = d["cross_validation_fraction"] if "cross_validation_fraction" in d else 0.25
        self.hyperparameters = d["hyperparameters"] if "hyperparameters" in d else None

def make_sdf_ref(fasta_ref, rtg, out):
    sp.call([rtg, 'format', '-o', out, fasta_ref])

def download_truth_set(name, reference, sample, out):
    ftp_vcf = known_truth_set_urls[name][reference][sample]['vcf']
    local_vcf = join(out, basename(ftp_vcf))
    urllib.request.urlretrieve(ftp_vcf, local_vcf)
    ftp_vcf_idx = known_truth_set_urls[name][reference][sample]['vcf_idx']
    local_vcf_idx = join(out, basename(ftp_vcf_idx))
    urllib.request.urlretrieve(ftp_vcf_idx, local_vcf_idx)
    ftp_bed = known_truth_set_urls[name][reference][sample]['bed']
    local_bed = join(out, basename(ftp_bed))
    urllib.request.urlretrieve(ftp_bed, local_bed)
    return local_vcf, local_vcf_idx, local_bed

def load_training_config(options):
    with open(options.config) as config:
        config = json.load(config)
        
        examples = config['examples'] if "examples" in config else config
        if type(examples) is not list:
            examples = [examples]
        examples = [TrainingData(example) for example in examples]

        assert all(example.reads is not None or example.octopus_vcf is not None for example in examples)
        check_exists([example.reference for example in examples])
        check_exists([example.regions for example in examples])
        check_exists_or_none([reads for example in examples for reads in example.reads])
        check_exists_or_none([example.octopus_vcf for example in examples])
        check_exists_or_none([example.config for example in examples])
        
        sdf_references, known_truths, given_truths = {}, {}, config["truths"] if "truths" in config else {}
        
        for label, details in given_truths.items():
            if not exists(details["vcf"]):
                raise ValueError(details["vcf"] + " does not exist")
            if not exists(details["bed"]):
                raise ValueError(details["bed"] + " does not exist")
        
        for example in examples:
            if example.sdf is None:
                if example.reference in sdf_references:
                    example.sdf = sdf_references[example.reference]
                else:
                    example.sdf = join(options.out, basename(example.reference).replace('fa', '').replace('fasta', '') + 'sdf')
                    make_sdf_ref(example.reference, options.rtg, example.sdf)
                sdf_references[example.reference] = example.sdf
            else:
                if not exists(example.sdf):
                    raise ValueError(example.sdf + " does not exist")
            
            if type(example.truth) is dict:
                if example.confident is None:
                    example.confident = {}
                for sample, truth_name in example.truth.items():
                    if truth_name in given_truths:
                        example.truth[sample], example.confident[sample] = given_truths[truth_name]['vcf'], given_truths[truth_name]['bed']
                    elif truth_name in known_truths:
                        example.truth[sample], example.confident[sample] = known_truths[example.truth]
                    else:
                        library, reference_version, truth_sample = truth_name.split('//')
                        vcf, _, bed = download_truth_set(library, reference_version, truth_sample, options.out)
                        known_truths[truth_name] = vcf, bed
                        example.truth[sample], example.confident[sample] = vcf, bed
            else:
                if not exists(example.truth):
                    if example.truth in given_truths:
                        example.truth, example.confident = given_truths[example.truth]['vcf'], given_truths[example.truth]['bed']
                    elif example.truth in known_truths:
                        example.truth, example.confident = known_truths[example.truth]
                    else:
                        name, reference_version, sample = example.truth.split('//')
                        vcf, _, bed = download_truth_set(name, reference_version, sample, options.out)
                        known_truths[example.truth] = vcf, bed
                        example.truth, example.confident = vcf, bed
                else:
                    if not exists(example.confident):
                        raise ValueError(example.confident + " does not exist")
                
        return examples, TrainingOptions(config['training']) if "training" in config else None

def get_reference_id(reference_filename):
    return basename(reference_filename).replace(".fasta", "")

def get_bam_id(bam_filenames):
    return "_".join(basename(fname).replace(".bam", "") for fname in bam_filenames)

def get_octopus_output_filename(reference_filename, bam_filenames, kind="germline"):
    return get_bam_id(bam_filenames) + "." + get_reference_id(reference_filename) + ".Octopus." + kind + ".vcf.gz"

def run_octopus(octopus, reference, reads, regions, threads, output,
                config=None, octopus_vcf=None, kind="germline"):
    octopus_cmd = [octopus, '-R', reference, '-I'] + reads + ['-t', regions,
                   '--ignore-unmapped-contigs', '--disable-call-filtering', '--annotations', 'forest',
                   '--threads', str(threads), '-o', output]
    if config is not None:
        octopus_cmd += ['--config', config]
    if octopus_vcf is not None:
        octopus_cmd += ['--filter-vcf', octopus_vcf]
    if kind == "somatic":
        octopus_cmd += ['--caller', 'cancer', '--somatics-only']
    sp.call(octopus_cmd)

def get_vcf_samples(vcf_filename):
    vcf = ps.VariantFile(vcf_filename)
    return vcf.header.samples

def run_rtg(rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, octopus_vcf_path, out_dir,
            bed_regions=None, sample=None, kind="germline"):
    cmd = [rtg, 'vcfeval', \
           '-t', rtg_ref_path, \
           '-b', truth_vcf_path, \
           '--evaluation-regions', confident_bed_path, \
           '-c', octopus_vcf_path, \
           '-o', out_dir]
    if bed_regions is not None:
        cmd += ['--bed-regions', bed_regions]
    if kind == "somatic":
        cmd += ['--squash-ploidy', '--sample', 'ALT']
    elif sample is not None:
        truth_samples = get_vcf_samples(truth_vcf_path)
        if len(truth_samples) > 1:
            raise Exception("More than one sample in truth " + truth_vcf_path)
        if sample == truth_samples[0]:
            cmd += ['--sample', sample]
        else:
            cmd += ['--sample', truth_samples[0] + "," + sample]
    sp.call(cmd)

def read_vcf_samples(vcf_fname):
    vcf = ps. VariantFile(vcf_fname)
    return vcf.header.samples

def get_annotation(field, rec, sample=None):
    if field == 'QUAL':
        return rec.qual
    elif field in rec.format:
        if sample is None:
            res = rec.samples[list(rec.samples)[0]][field]
        else:
            res = rec.samples[sample][field]
        if type(res) == tuple:
            res = res[0]
        return res
    else:
        res = rec.info[field]
        if type(res) == tuple:
            res = res[0]
        return res

def is_somatic_sample(rec, sample):
    return bool(int(get_annotation('SOMATIC', rec, sample)))

def is_somatic_record(rec):
    return any(is_somatic_sample(rec, sample) for sample in list(rec.samples))

def index_vcf(vcf_filename):
    sp.call(['tabix', '-f', vcf_filename])

def filter_somatic(vcf_filename):
    in_vcf = ps. VariantFile(vcf_filename)
    tmp_vcf_filename = vcf_filename.replace('.vcf', '.tmp.vcf')
    out_vcf = ps.VariantFile(tmp_vcf_filename, 'wz', header=in_vcf.header)
    num_skipped_records = 0
    for rec in in_vcf:
        if is_somatic_record(rec):
            try:
                out_vcf.write(rec)
            except OSError:
                num_skipped_records += 1
    if num_skipped_records > 0:
        print("Skipped " + str(num_skipped_records) + " bad records")
    in_vcf.close()
    out_vcf.close()
    shutil.move(tmp_vcf_filename, vcf_filename)
    index_vcf(vcf_filename)

def read_octopus_header_info(vcf_filename):
    vcf = ps.VariantFile(vcf_filename)
    for record in vcf.header.records:
        if record.key == "octopus":
            return dict(record)
    return None

def read_normal_samples(vcf_filename):
    options = read_octopus_header_info(vcf_filename)['options'].split(' ')
    result = []
    for token in options[options.index('--normal-sample') + 1:]:
        if token.startswith('--'):
            break
        else:
            result.append(token)
    return result

def is_normal_sample(sample, vcf_filename):
    return sample in read_normal_samples(vcf_filename)

def eval_octopus(octopus, rtg, example, out_dir, threads, kind="germline"):
    if example.reads is not None:
        octopus_vcf = join(out_dir, get_octopus_output_filename(example.reference, example.reads, kind=kind))
        run_octopus(octopus, example.reference, example.reads, example.regions, threads, octopus_vcf,
                    config=example.config, octopus_vcf=example.octopus_vcf, kind=kind)
    else:
        assert example.octopus_vcf is not None
        octopus_vcf = example.octopus_vcf
    if kind == "somatic":
        # Hack as '--somatics-only' option is currently ignored when in training mode
        filter_somatic(octopus_vcf)
    samples = read_vcf_samples(octopus_vcf)
    result = []
    if len(samples) == 1:
        rtf_eval_dir = join(out_dir, basename(octopus_vcf).replace(".vcf.gz", ".") + '.eval')
        run_rtg(rtg, example.sdf, example.truth, example.confident, octopus_vcf, rtf_eval_dir, bed_regions=example.regions, kind=kind)
        result.append(rtf_eval_dir)
    elif kind == "somatic":
        for sample in samples:
            if not is_normal_sample(sample, octopus_vcf):
                rtf_eval_dir = join(out_dir, basename(octopus_vcf).replace(".vcf.gz", ".") + sample + '.eval')
                run_rtg(rtg, example.sdf, example.truth, example.confident, octopus_vcf, rtf_eval_dir, bed_regions=example.regions, sample=sample, kind=kind)
                result.append(rtf_eval_dir)
    else:
        for sample in samples:
            if sample in example.truth:
                rtf_eval_dir = join(out_dir, basename(octopus_vcf).replace(".vcf.gz", ".") + sample + '.eval')
                run_rtg(rtg, example.sdf, example.truth[sample], example.confident[sample], octopus_vcf, rtf_eval_dir, bed_regions=example.regions, sample=sample, kind=kind)
                result.append(rtf_eval_dir)
    return result

def subset(vcf_in_path, vcf_out_path, bed_regions):
    sp.call(['bcftools', 'view', '-R', bed_regions, '-O', 'z', '-o', vcf_out_path, vcf_in_path])

def is_missing(x):
    return x == '.' or np.isnan(float(x))

def annotation_to_string(x, missing_value):
    return str(missing_value) if is_missing(x) else str(x)

def make_ranger_data(octopus_vcf_path, out_path, classifcation, measures, missing_value=-1, fraction=1):
    vcf = ps. VariantFile(octopus_vcf_path)
    with open(out_path, 'w') as ranger_data:
        datawriter = csv.writer(ranger_data, delimiter=' ')
        for rec in vcf:
            if fraction >= 1 or random.random() <= fraction:
                row = [annotation_to_string(get_annotation(measure, rec), missing_value) for measure in measures]
                row.append(str(int(classifcation)))
                datawriter.writerow(row)

def concat(filenames, outpath):
    with open(outpath, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def shuffle(fname):
    lines = open(fname).readlines()
    random.shuffle(lines)
    open(fname, 'w').writelines(lines)

def add_header(fname, header):
    lines = open(fname).readlines()
    with open(fname, 'w') as f:
        f.write(header + '\n')
        f.writelines(lines)

def partition_data(data_fname, validation_fraction, training_fname, validation_fname):
    with open(data_fname) as data_file:
        header = data_file.readline()
        with open(training_fname, 'w') as training_file, open(validation_fname, 'w') as validation_file:
            training_file.write(header)
            validation_file.write(header)
            for example in data_file:
                if random.random() < validation_fraction:
                    validation_file.write(example)
                else:
                    training_file.write(example)

def run_ranger_training(ranger, data_path, hyperparameters, threads, out, seed=None):
    cmd = [ranger, '--file', data_path, '--depvarname', 'TP', '--probability',
          '--nthreads', str(threads), '--outprefix', out, '--write', '--impmeasure', '1', '--verbose']
    if 'trees' in hyperparameters:
        cmd += ['--ntree', str(hyperparameters["trees"])]
    if 'targetpartitionsize' in hyperparameters:
        cmd += ['--targetpartitionsize', str(hyperparameters["min_node_size"])]
    if 'maxdepth' in hyperparameters:
        cmd += ['--maxdepth', str(hyperparameters["maxdepth"])]
    if seed is not None:
        cmd += ['--seed', str(seed)]
    sp.call(cmd)

def run_ranger_prediction(ranger, forest, data_path, threads, out):
    cmd = [ranger, '--file', data_path, '--predict', forest,
          '--nthreads', str(threads), '--outprefix', out, '--verbose']
    sp.call(cmd)

def read_predictions(prediction_fname):
    with open(prediction_fname) as prediction_file:
        next(prediction_file)
        next(prediction_file)
        next(prediction_file)
        return np.array([float(line.strip().split()[0]) for line in prediction_file])

def read_truth_labels(truth_fname):
    with open(truth_fname) as truth_file:
        true_label_index = truth_file.readline().strip().split().index('TP')
        return np.array([int(line.strip().split()[true_label_index]) for line in truth_file])

def select_training_hypterparameters(master_data_fname, options):
    if options is None or options.hyperparameters is None:
        return {"trees": 500, "min_node_size": 10}
    elif len(options.hyperparameters) == 1:
        return options.hyperparameters[0]
    else:
        # Use cross-validation to select hypterparameters
        training_data_fname, validation_data_fname = master_data_fname.replace('.dat', '.train.data'), master_data_fname.replace('.dat', '.validate.data')
        partition_data(master_data_fname, options.cross_validation_fraction, training_data_fname, validation_data_fname)
        optimal_params, min_loss = None, None
        cross_validation_temp_dir = join(options.out, 'cross_validation')
        makedirs(cross_validation_temp_dir)
        cross_validation_prefix = join(cross_validation_temp_dir, options.prefix)
        prediction_fname = cross_validation_prefix + '.prediction'
        forest_fname = cross_validation_prefix + '.forest'
        for params in options.hyperparameters:
            trees, min_node_size = params["trees"], params["min_node_size"]
            print('Training cross validation forest with trees =', trees, 'min_node_size =', min_node_size)
            run_ranger_training(options.ranger, training_data_fname, trees, min_node_size, options.threads, cross_validation_prefix, seed=10)
            run_ranger_prediction(options.ranger, forest_fname, validation_data_fname, options.threads, cross_validation_prefix)
            truth_labels, predictions = read_truth_labels(validation_data_fname), read_predictions(prediction_fname)
            loss = log_loss(truth_labels, predictions)
            print('Binary cross entropy =', loss)
            if min_loss is None or loss < min_loss:
                min_loss, optimal_params = loss, params
        shutil.rmtree(cross_validation_temp_dir)
        return optimal_params

def main(options):
    examples, training_params = load_training_config(options)
    if not exists(options.out):
        makedirs(options.out)
    data_files, tmp_files = [], []
    for example in examples:
        rtg_eval_dirs = eval_octopus(options.octopus, options.rtg, example, options.out, options.threads, kind=options.kind)
        for rtg_eval in rtg_eval_dirs:
            tp_vcf_path = join(rtg_eval, "tp.vcf.gz")
            tp_train_vcf_path = tp_vcf_path.replace("tp.vcf", "tp.train.vcf")
            subset(tp_vcf_path, tp_train_vcf_path, example.regions)
            tp_data_path = tp_train_vcf_path.replace(".vcf.gz", ".dat")
            make_ranger_data(tp_train_vcf_path, tp_data_path, True, default_measures, options.missing_value, fraction=example.tp)
            data_files.append(tp_data_path)
            fp_vcf_path = join(rtg_eval, "fp.vcf.gz")
            fp_train_vcf_path = fp_vcf_path.replace("fp.vcf", "fp.train.vcf")
            subset(fp_vcf_path, fp_train_vcf_path, example.regions)
            fp_data_path = fp_train_vcf_path.replace(".vcf.gz", ".dat")
            make_ranger_data(fp_train_vcf_path, fp_data_path, False, default_measures, options.missing_value, fraction=example.fp)
            data_files.append(fp_data_path)
            tmp_files += [tp_train_vcf_path, fp_train_vcf_path]
    master_data_file = join(options.out, options.prefix + ".dat")
    concat(data_files, master_data_file)
    for file in tmp_files + data_files:
        remove(file)
    shuffle(master_data_file)
    ranger_header = ' '.join(default_measures + ['TP'])
    add_header(master_data_file, ranger_header)
    ranger_out_prefix = join(options.out, options.prefix)
    hyperparameters = select_training_hypterparameters(master_data_file, training_params)
    run_ranger_training(options.ranger, master_data_file, hyperparameters, options.threads, ranger_out_prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',
                        required=True,
                        type=str,
                        help='Training data config in json format')
    parser.add_argument('--octopus',
                        type=str,
                        default=default_octopus_bin,
                        help='Octopus binary')
    parser.add_argument('--rtg', 
                        type=str,
                        default='rtg',
                        help='RTG Tools binary')
    parser.add_argument('--ranger', 
                        type=str,
                        default='ranger',
                        help='Ranger binary')
    parser.add_argument('-o', '--out',
                        type=str,
                        required=True,
                        help='Output directory')
    parser.add_argument('--prefix',
                        type=str,
                        default='octopus',
                        help='Output files prefix')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of threads for octopus')
    parser.add_argument('--missing_value',
                        type=float,
                        default=-1,
                        help='Value for missing measures')
    parser.add_argument('--kind',
                        type=str,
                        default='germline',
                        help='Kind of random forest to train [germline, somatic]')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
