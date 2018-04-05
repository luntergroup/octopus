#!/usr/bin/env python3

import argparse
from os import makedirs, remove
from os.path import join, basename, exists
from subprocess import call
import csv
from pysam import VariantFile
import random

def run_octopus(octopus, ref_path, bam_path, regions_bed, measures, threads, out_path):
    call([octopus, '-R', ref_path, '-I', bam_path, '-t', regions_bed, '-o', out_path, '--threads', str(threads), '--legacy', '--csr-train'] + measures)

def get_reference_id(ref_path):
    return basename(ref_path).replace(".fasta", "")

def get_bam_id(bam_path):
    return basename(bam_path).replace(".bam", "")

def call_variants(octopus, ref_path, bam_path, regions_bed, measures, threads, out_dir):
    ref_id = get_reference_id(ref_path)
    bam_id = get_bam_id(bam_path)
    out_vcf = join(out_dir, "octopus." + bam_id + "." + ref_id + ".vcf.gz")
    run_octopus(octopus, ref_path, bam_path, regions_bed, measures, threads, out_vcf)
    legacy_vcf = out_vcf.replace(".vcf.gz", ".legacy.vcf.gz")
    return legacy_vcf

def run_rtg(rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, octopus_vcf_path, out_dir):
    call([rtg, 'vcfeval', '-b', truth_vcf_path, '-t', rtg_ref_path, '--evaluation-regions', confident_bed_path, '--ref-overlap', '-c', octopus_vcf_path, '-o', out_dir])

def eval_octopus(octopus, ref_path, bam_path, regions_bed, measures, threads, rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, out_dir):
    octopus_vcf = call_variants(octopus, ref_path, bam_path, regions_bed, measures, threads, out_dir)
    rtf_eval_dir = join(out_dir, basename(octopus_vcf).replace(".legacy.vcf.gz", ".eval"))
    run_rtg(rtg, rtg_ref_path, truth_vcf_path, confident_bed_path, octopus_vcf, rtf_eval_dir)
    return rtf_eval_dir

def is_missing(x):
    return x == '.'

def to_str(x, missing_value=0):
    if is_missing(x):
        return str(missing_value)
    else:
        return str(x)

def get_field(field, rec, missing_value=0):
    if field == 'QUAL':
        return to_str(rec.qual, missing_value)
    else:
        val = rec.info[field]
        if type(val) == tuple:
            val = val[0]
        return to_str(val, missing_value)

def make_ranger_data(octopus_vcf_path, measures, is_tp, out, missing_value=0):
    vcf = VariantFile(octopus_vcf_path)
    with open(out, 'w') as ranger_dat:
        datwriter = csv.writer(ranger_dat, delimiter=' ')
        for rec in vcf:
            row = [get_field(measure, rec, missing_value) for measure in measures]
            if is_tp:
                row.append('1')
            else:
                row.append('0')
            datwriter.writerow(row)

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

def run_ranger_training(ranger, data_path, threads, out):
    call([ranger, '--file', data_path, '--depvarname', 'TP', '--probability', '--nthreads', str(threads), '--outprefix', out, '--write', out + '.forest'])

def main(options):
    if not exists(options.out):
        makedirs(options.out)
    rtg_eval_dirs = []
    for bam_path in options.reads:
        rtg_eval_dirs.append(eval_octopus(options.octopus, options.reference, bam_path, options.regions, options.measures, options.threads,
                                      options.rtg, options.sdf, options.truth, options.confident, options.out))
    data_paths = []
    for rtg_eval in rtg_eval_dirs:
        tp_vcf_path = join(rtg_eval, "tp.vcf.gz")
        tp_data_path = tp_vcf_path.replace(".vcf.gz", ".dat")
        make_ranger_data(tp_vcf_path, options.measures, True, tp_data_path)
        data_paths.append(tp_data_path)
        fp_vcf_path = join(rtg_eval, "fp.vcf.gz")
        fp_data_path = fp_vcf_path.replace(".vcf.gz", ".dat")
        make_ranger_data(fp_vcf_path, options.measures, False, fp_data_path)
        data_paths.append(fp_data_path)
    master_data_path = join(options.out, "ranger_train_master.dat")
    concat(data_paths, master_data_path)
    for path in data_paths:
        remove(path)
    shuffle(master_data_path)
    ranger_header = ' '.join(options.measures + ['TP'])
    add_header(master_data_path, ranger_header)
    ranger_out_prefix = join(options.out, "ranger_octopus")
    run_ranger_training(options.ranger, master_data_path, options.threads, ranger_out_prefix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-R', '--reference',
                        type=str,
                        required=True,
                        help='Reference to use for calling')
    parser.add_argument('-I', '--reads',
                        nargs='+',
                        type=str,
                        required=True,
                        help='Input BAM files')
    parser.add_argument('--regions',
                        type=str,
                        required=True,
                        help='BED files containing regions to call')
    parser.add_argument('--measures',
                        type=str,
                        nargs='+',
                        required=True,
                        help='Measures to use for training')
    parser.add_argument('--truth',
                        type=str,
                        required=True,
                        help='Truth VCF file')
    parser.add_argument('--confident',
                        type=str,
                        required=True,
                        help='BED files containing high confidence truth regions')
    parser.add_argument('--octopus',
                        type=str,
                        required=True,
                        help='Octopus binary')
    parser.add_argument('--rtg', 
                        type=str,
                        required=True,
                        help='RTG Tools binary')
    parser.add_argument('--sdf', 
                        type=str,
                        required=True,
                        help='RTG Tools SDF reference index')
    parser.add_argument('--ranger', 
                        type=str,
                        required=True,
                        help='Ranger binary')
    parser.add_argument('-o', '--out',
                        type=str,
                        help='Output directory')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of threads for octopus')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
