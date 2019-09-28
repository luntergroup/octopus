#!/usr/bin/env python3

import argparse
from os import makedirs, remove, pardir
from os.path import join, basename, exists, dirname, abspath
from subprocess import call
import csv
from pysam import VariantFile
import random
import numpy as np
import shutil
from sklearn.metrics import log_loss
import json
import urllib.request

script_dir = dirname(abspath(__file__))
default_octopus_bin = join(abspath(join(script_dir, pardir)), 'bin/octopus')

default_measures = "AC AD ADP AF ARF BQ CC CRF DP ER ERS FRF GC GQ GQD NC MC MF MP MRC MQ MQ0 MQD PP PPD QD QUAL REFCALL REB RSB RTB SB SD SF SHC SMQ SOMATIC STR_LENGTH STR_PERIOD VL".split()

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

def check_exists(paths):
    for path in paths:
        if not exists(path):
            raise ValueError(path + " does not exist")

def check_exists_or_none(paths):
    for path in paths:
        if path is not None and not exists(path):
            raise ValueError(path + " does not exist")

class TrainingData:
    def __init__(self, data):
        self.reference = data["reference"] if "reference" in data else None
        self.sdf = data["sdf"] if "sdf" in data else None
        self.reads = data["reads"] if "reads" in data else None
        self.regions = data["regions"] if "regions" in data else None
        self.truth = data["truth"] if "truth" in data else None
        self.confident = data["confident"] if "confident" in data else None
        self.tp = data["tp"] if "tp" in data else 1
        self.fp = data["fp"] if "tp" in data else 1
        self.config = data["config"] if "config" in data else None

def read_training_data(training_json_filename):
    with open(training_json_filename) as training_json:
        return [TrainingData(data) for data in json.load(training_json)]

def make_sdf_ref(fasta_ref, rtg, out):
    call([rtg, 'format', '-o', out, fasta_ref])

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

def setup(training_configs, rtg, out):
    check_exists([config.reference for config in training_configs])
    check_exists([config.reads for config in training_configs])
    check_exists([config.regions for config in training_configs])
    check_exists_or_none([config.config for config in training_configs])
    result = []
    sdf_references, known_truths = {}, {}
    for data in training_configs:
        if data.sdf is None:
            if data.reference in sdf_references:
                data.sdf = sdf_references[data.reference]
            else:
                data.sdf = join(out, basename(data.reference).replace('fa', '').replace('fasta', '') + 'sdf')
                make_sdf_ref(data.reference, rtg, data.sdf)
            sdf_references[data.reference] = data.sdf
        else:
            if not exists(data.sdf):
                raise ValueError(data.sdf + " does not exist")
        if not exists(data.truth):
            if data.truth in known_truths:
                data.truth, data.confident = known_truths[data.truth]
            else:
                try:
                    name, reference_version, sample = data.truth.split('.')
                    vcf, idx, bed = download_truth_set(name, reference_version, sample, out)
                    known_truths[data.truth] = vcf, bed
                    data.truth, data.confident = vcf, bed
                except:
                    print('Invalid truth set format')
        else:
            if not exists(data.confident):
                raise ValueError(data.confident + " does not exist")
        result.append(data)
    return result

def load_training_data(options):
    return setup(read_training_data(options.config), options.rtg, options.out)

def run_octopus(octopus, ref_path, bam_path, regions_bed, threads, out_path, config=None):
    octopus_cmd = [octopus, '-R', ref_path, '-I', bam_path, '-t', regions_bed,
                   '--ignore-unmapped-contigs', '--disable-call-filtering', '--annotations', 'forest',
                   '--threads', str(threads), '-o', out_path]
    if config is not None:
        octopus_cmd += ['--config', config]
    call(octopus_cmd)

def get_reference_id(ref_path):
    return basename(ref_path).replace(".fasta", "")

def get_bam_id(bam_path):
    return basename(bam_path).replace(".bam", "")

def call_variants(octopus, ref_path, bam_path, regions_bed, threads, out_dir, config=None):
    ref_id = get_reference_id(ref_path)
    bam_id = get_bam_id(bam_path)
    out_vcf = join(out_dir, "octopus." + bam_id + "." + ref_id + ".vcf.gz")
    run_octopus(octopus, ref_path, bam_path, regions_bed, threads, out_vcf, config=config)
    return out_vcf

def run_rtg(rtg, rtg_ref_path, truth_vcf_path, bed_regions, confident_bed_path, octopus_vcf_path, out_dir):
    call([rtg, 'vcfeval', '-b', truth_vcf_path, '-t', rtg_ref_path,
         '--bed-regions', bed_regions,
         '--evaluation-regions', confident_bed_path,
          '--ref-overlap', '-c', octopus_vcf_path, '-o', out_dir])

def eval_octopus(octopus, rtg, training_data, out_dir, threads):
    octopus_vcf = call_variants(octopus, training_data.reference, training_data.reads, training_data.regions, threads, out_dir, config=training_data.config)
    rtf_eval_dir = join(out_dir, basename(octopus_vcf).replace(".vcf.gz", ".eval"))
    run_rtg(rtg, training_data.sdf, training_data.truth, training_data.regions, training_data.confident, octopus_vcf, rtf_eval_dir)
    return rtf_eval_dir

def subset(vcf_in_path, vcf_out_path, bed_regions):
    call(['bcftools', 'view', '-R', bed_regions, '-O', 'z', '-o', vcf_out_path, vcf_in_path])

def get_annotation(field, rec):
    if field == 'QUAL':
        return rec.qual
    elif field in rec.format:
        res = rec.samples[list(rec.samples)[0]][field]
        if type(res) == tuple:
            res = res[0]
        return res
    else:
        res = rec.info[field]
        if type(res) == tuple:
            res = res[0]
        return res

def is_missing(x):
    return x == '.' or np.isnan(float(x))

def annotation_to_string(x, missing_value):
    return str(missing_value) if is_missing(x) else str(x)

def make_ranger_data(octopus_vcf_path, out_path, classifcation, measures, missing_value=-1, fraction=1):
    vcf = VariantFile(octopus_vcf_path)
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

def run_ranger_training(ranger, data_path, n_trees, min_node_size, threads, out, seed=None):
    cmd = [ranger, '--file', data_path, '--depvarname', 'TP', '--probability',
          '--ntree', str(n_trees), '--targetpartitionsize', str(min_node_size),
          '--nthreads', str(threads), '--outprefix', out, '--write', '--verbose']
    if seed is not None:
        cmd += ['--seed', str(seed)]
    call(cmd)

def run_ranger_prediction(ranger, forest, data_path, threads, out):
    call([ranger, '--file', data_path, '--predict', forest,
          '--nthreads', str(threads), '--outprefix', out, '--verbose'])

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
    if type(options.trees) is int and type(options.min_node_size) is int:
        return options.trees, options.min_node_size
    else:
        # Use cross-validation to select hypterparameters
        training_data_fname, validation_data_fname = master_data_fname.replace('.dat', '.train.data'), master_data_fname.replace('.dat', '.validate.data')
        partition_data(master_data_fname, options.cross_validation_fraction, training_data_fname, validation_data_fname)
        optimal_params = None, None
        min_loss = None
        cross_validation_temp_dir = join(options.out, 'cross_validation')
        makedirs(cross_validation_temp_dir)
        cross_validation_prefix = join(cross_validation_temp_dir, options.prefix)
        prediction_fname = cross_validation_prefix + '.prediction'
        forest_fname = cross_validation_prefix + '.forest'
        for trees in ([options.trees] if type(options.trees) is int else options.trees):
            for min_node_size in ([options.min_node_size] if type(options.min_node_size) is int else options.min_node_size):
                print('Training cross validation forest with trees =', trees, 'min_node_size =', min_node_size)
                run_ranger_training(options.ranger, training_data_fname, trees, min_node_size, options.threads, cross_validation_prefix, seed=10)
                run_ranger_prediction(options.ranger, forest_fname, validation_data_fname, options.threads, cross_validation_prefix)
                truth_labels, predictions = read_truth_labels(validation_data_fname), read_predictions(prediction_fname)
                loss = log_loss(truth_labels, predictions)
                print('Binary cross entropy =', loss)
                if min_loss is None or loss < min_loss:
                    min_loss = loss
                    optimal_params = trees, min_node_size
        shutil.rmtree(cross_validation_temp_dir)
        return optimal_params

def main(options):
    if not exists(options.out):
        makedirs(options.out)
    data_files, tmp_files = [], []
    for config in load_training_data(options):
        rtg_eval = eval_octopus(options.octopus, options.rtg, config, options.out, options.threads)
        tp_vcf_path = join(rtg_eval, "tp.vcf.gz")
        tp_train_vcf_path = tp_vcf_path.replace("tp.vcf", "tp.train.vcf")
        subset(tp_vcf_path, tp_train_vcf_path, config.regions)
        tp_data_path = tp_train_vcf_path.replace(".vcf.gz", ".dat")
        make_ranger_data(tp_train_vcf_path, tp_data_path, True, default_measures, options.missing_value, fraction=config.tp)
        data_files.append(tp_data_path)
        fp_vcf_path = join(rtg_eval, "fp.vcf.gz")
        fp_train_vcf_path = fp_vcf_path.replace("fp.vcf", "fp.train.vcf")
        subset(fp_vcf_path, fp_train_vcf_path, config.regions)
        fp_data_path = fp_train_vcf_path.replace(".vcf.gz", ".dat")
        make_ranger_data(fp_train_vcf_path, fp_data_path, False, default_measures, options.missing_value, fraction=config.fp)
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
    num_trees, min_node_size = select_training_hypterparameters(master_data_file, options)
    run_ranger_training(options.ranger, master_data_file, num_trees, min_node_size, options.threads, ranger_out_prefix)
    try:
        remove(ranger_out_prefix + ".confusion")
    except FileNotFoundError:
        print('Ranger did not complete, perhaps it was killed due to insufficient memory')

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
                        required=True,
                        help='RTG Tools binary')
    parser.add_argument('--ranger', 
                        type=str,
                        required=True,
                        help='Ranger binary')
    parser.add_argument('--trees',
                        nargs='+',
                        type=int,
                        default=300,
                        help='Number of trees to use in the random forest')
    parser.add_argument('--min_node_size',
                        nargs='+',
                        type=int,
                        default=20,
                        help='Node size to stop growing trees, implicitly limiting tree depth')
    parser.add_argument('--cross_validation_fraction',
                        type=float,
                        default=0.25,
                        help='Fraction of data points to hold back for cross validation')
    parser.add_argument('-o', '--out',
                        type=str,
                        required=True,
                        help='Output directory')
    parser.add_argument('--prefix',
                        type=str,
                        default='octopus_germline',
                        help='Output files prefix')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of threads for octopus')
    parser.add_argument('--missing_value',
                        type=float,
                        default=-1,
                        help='Value for missing measures')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
