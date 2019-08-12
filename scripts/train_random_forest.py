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

script_dir = dirname(abspath(__file__))
default_octopus_bin = join(abspath(join(script_dir, pardir)), 'bin/octopus')

default_measures = "AC AD ADP AF ARF BQ CC CRF DP FRF GC GQ GQD NC MC MF MP MRC MQ MQ0 MQD PP PPD QD QUAL REFCALL REB RSB RTB SB SD SF SHC SMQ SOMATIC STR_LENGTH STR_PERIOD VL".split()

def make_sdf_ref(fasta_ref, rtg, out):
    call([rtg, 'format', '-o', out, fasta_ref])

def run_octopus(octopus, ref_path, bam_path, regions_bed, threads, out_path, config=None):
    octopus_cmd = [octopus, '-R', ref_path, '-I', bam_path, '-t', regions_bed,
                   '-f', 'off', '--annotations', 'forest',
                   '--legacy', '--threads', str(threads), '-o', out_path]
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
    legacy_vcf = out_vcf.replace(".vcf.gz", ".legacy.vcf.gz")
    return legacy_vcf

def run_rtg(rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, octopus_vcf_path, out_dir):
    call([rtg, 'vcfeval', '-b', truth_vcf_path, '-t', rtg_ref_path, '--evaluation-regions', confident_bed_path,
          '--ref-overlap', '-c', octopus_vcf_path, '-o', out_dir])

def eval_octopus(octopus, ref_path, bam_path, regions_bed, threads,
                 rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, out_dir,
                 config=None):
    octopus_vcf = call_variants(octopus, ref_path, bam_path, regions_bed, threads, out_dir, config=config)
    rtf_eval_dir = join(out_dir, basename(octopus_vcf).replace(".legacy.vcf.gz", ".eval"))
    run_rtg(rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, octopus_vcf, rtf_eval_dir)
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

def make_ranger_data(octopus_vcf_path, out_path, classifcation, measures, missing_value=-1):
    vcf = VariantFile(octopus_vcf_path)
    with open(out_path, 'w') as ranger_data:
        datawriter = csv.writer(ranger_data, delimiter=' ')
        for rec in vcf:
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

def run_ranger_training(ranger, data_path, n_trees, min_node_size, threads, out):
    call([ranger, '--file', data_path, '--depvarname', 'TP', '--probability',
          '--ntree', str(n_trees), '--targetpartitionsize', str(min_node_size),
          '--nthreads', str(threads), '--outprefix', out, '--write', '--verbose'])

def main(options):
    if not exists(options.out):
        makedirs(options.out)
    tmp_files, tmp_dirs, rtg_eval_dirs = [], [], []
    truth_sets, confident_sets = options.truth, options.confident
    fasta_refs, rtg_sdf_refs = options.reference, []
    call_region_beds = options.regions
    if len(call_region_beds) == 1:
        call_region_beds = len(options.reads) * [call_region_beds[0]]
    for ref in fasta_refs:
        rtg_sdf_ref = basename(ref) + '.tmp.sdf'
        make_sdf_ref(ref, options.rtg, rtg_sdf_ref)
        rtg_sdf_refs.append(rtg_sdf_ref)
        tmp_dirs.append(rtg_sdf_ref)
    if len(options.reference) == 1 and len(options.reads) > 1:
        fasta_refs = len(options.reads) * [fasta_refs[0]]
        rtg_sdf_refs = len(options.reads) * [rtg_sdf_refs[0]]
    if len(truth_sets) == 1 and len(options.reads) > 1:
        truth_sets, confident_sets = len(options.reads) * [truth_sets[0]], len(options.reads) * [confident_sets[0]]
    configs = options.config
    if configs is None or len(configs) == 1 and len(options.reads) > 1:
        configs = len(options.reads) * [configs]
    for fasta_ref, bam_path, regions_bed, config, rtg_ref, truth, confident in \
            zip(fasta_refs, options.reads, call_region_beds, configs, rtg_sdf_refs, truth_sets, confident_sets):
        rtg_eval_dirs.append(eval_octopus(options.octopus, fasta_ref, bam_path, regions_bed,
                                          options.threads, options.rtg, rtg_ref, truth, confident,
                                          options.out, config=config))
    data_files = []
    for rtg_eval, regions_bed in zip(rtg_eval_dirs, call_region_beds):
        tp_vcf_path = join(rtg_eval, "tp.vcf.gz")
        tp_train_vcf_path = tp_vcf_path.replace("tp.vcf", "tp.train.vcf")
        subset(tp_vcf_path, tp_train_vcf_path, regions_bed)
        tp_data_path = tp_train_vcf_path.replace(".vcf.gz", ".dat")
        make_ranger_data(tp_train_vcf_path, tp_data_path, True, default_measures, options.missing_value)
        data_files.append(tp_data_path)
        fp_vcf_path = join(rtg_eval, "fp.vcf.gz")
        fp_train_vcf_path = fp_vcf_path.replace("fp.vcf", "fp.train.vcf")
        subset(fp_vcf_path, fp_train_vcf_path, regions_bed)
        fp_data_path = fp_train_vcf_path.replace(".vcf.gz", ".dat")
        make_ranger_data(fp_train_vcf_path, fp_data_path, False, default_measures, options.missing_value)
        data_files.append(fp_data_path)
        tmp_files += [tp_train_vcf_path, fp_train_vcf_path]
    master_data_file = join(options.out, options.prefix + ".dat")
    concat(data_files, master_data_file)
    for path in data_files:
        remove(path)
    for file in tmp_files:
        remove(file)
    for dir in tmp_dirs:
        shutil.rmtree(dir)
    shuffle(master_data_file)
    ranger_header = ' '.join(default_measures + ['TP'])
    add_header(master_data_file, ranger_header)
    ranger_out_prefix = join(options.out, options.prefix)
    run_ranger_training(options.ranger, master_data_file, options.trees, options.min_node_size, options.threads, ranger_out_prefix)
    try:
        remove(ranger_out_prefix + ".confusion")
    except FileNotFoundError:
        print('Ranger did not complete, perhaps it was killed due to insufficient memory')

def check_options(options):
    if len(options.truth) > 1:
        if len(options.truth) != len(options.reads):
            print("Must specify a --truth set for each --read set if using different truth sets")
            return False
        if len(options.confident) != len(options.truth):
            print("Must specify a --confident set for each --truth set")
            return False
    if options.config is not None and len(options.config) > 1:
        if len(options.config) != len(options.reads):
            print("Must specify a --config file for each --read set if using different configs")
            return False
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-R', '--reference',
                        nargs='+',
                        type=str,
                        required=True,
                        help='Reference to use for calling')
    parser.add_argument('-I', '--reads',
                        nargs='+',
                        type=str,
                        required=True,
                        help='Input BAM files')
    parser.add_argument('-T', '--regions',
                        nargs='+',
                        type=str,
                        required=True,
                        help='BED files containing regions to call')
    parser.add_argument('--config',
                        nargs='+',
                        type=str,
                        help='Octopus config files for each read set, or one for all')
    parser.add_argument('--truth',
                        nargs='+',
                        type=str,
                        required=True,
                        help='Truth VCF file')
    parser.add_argument('--confident',
                        nargs='+',
                        type=str,
                        required=True,
                        help='BED files containing high confidence truth regions')
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
                        type=int,
                        default=300,
                        help='Number of trees to use in the random forest')
    parser.add_argument('--min_node_size',
                        type=int,
                        default=20,
                        help='Node size to stop growing trees, implicitly limiting tree depth')
    parser.add_argument('-o', '--out',
                        type=str,
                        required=True,
                        help='Output directory')
    parser.add_argument('--prefix',
                        type=str,
                        default='ranger_octopus',
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
    if not check_options(parsed):
        exit(1)
    main(parsed)
