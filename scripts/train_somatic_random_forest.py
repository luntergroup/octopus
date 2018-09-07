#!/usr/bin/env python3

import argparse
from os import makedirs, remove
from os.path import join, basename, exists
import shutil
from subprocess import call
import csv
from pysam import VariantFile
import random
import numpy as np

def get_sample_names(truth_vcf_name):
    truth_vcf = VariantFile(truth_vcf_name)
    return [sample for sample in truth_vcf.header.samples]

def index(vcf_path):
    call(['tabix', vcf_path])
    return vcf_path + '.tbi'

def remove_vcf(vcf_path):
    remove(vcf_path)
    if exists(vcf_path + '.tbi'):
        remove(vcf_path + '.tbi')

def is_homozygous(gt):
    return all(a == 1 for a in gt)

def to_float(val):
    if val == '.':
        return np.nan
    else:
        try:
            return float(val)
        except ValueError:
            return val

def get_value(rec, feature):
    if feature == "POS":
        return rec.pos
    elif feature == "REF":
        return rec.ref
    elif feature == "ALT":
        return list(rec.alts)
    elif feature == "QUAL":
        return rec.qual
    elif feature == 'AS':
        return max([len(allele) for allele in rec.alleles])
    elif feature in rec.info:
        val = rec.info[feature]
        if type(val) == tuple:
            return tuple([to_float(v) for v in val])
        else:
            return to_float(val)
    elif feature in rec.samples[0]:
        if feature == 'GT':
            return int(is_homozygous(rec.samples[0][feature]))
        else:
            return int(rec.samples[0][feature])
    else:
        return np.nan

def is_somatic(rec):
    return any(get_value(rec, 'SOMATIC'))

def filter_somatic(in_vcf_path, out_vcf_path):
    in_vcf = VariantFile(in_vcf_path)
    out_vcf = VariantFile(out_vcf_path, 'w', header=in_vcf.header)
    num_skipped_records = 0
    for rec in in_vcf:
        if is_somatic(rec):
            try:
                out_vcf.write(rec)
            except OSError:
                num_skipped_records += 1
    print("Skipped " + str(num_skipped_records) + " bad records")
    in_vcf.close()
    out_vcf.close()

def prepare_octopus_vcf_for_rtg(octopus_vcf, tumour_sample, out_vcf_name):
    """"
    Octopus reports non-diploid genotypes for somatic variants.
    """
    in_vcf = VariantFile(octopus_vcf)
    out_vcf = VariantFile(out_vcf_name, 'w', header=in_vcf.header)
    n_failed = 0
    for record in in_vcf:
        old_gt = record.samples[tumour_sample]['GT']
        assert(len(old_gt) > 1)
        somatic_allele = next(a for a in reversed(list(old_gt)) if a is not None and a > 0)
        record.samples[tumour_sample]['GT'] = (old_gt[0], somatic_allele)
        try:
            out_vcf.write(record)
        except OSError:
            n_failed += 1
    out_vcf.close()
    index(out_vcf_name)

def run_rtg_vcfeval(rtg_bin, rtg_ref, caller_vcf, truth_vcf, sample, out_dir, regions_bed=None):
    rtg_cmd = [rtg_bin, 'vcfeval', '-t', rtg_ref, '-b', truth_vcf, '-c', caller_vcf, '-o', out_dir,
               '--squash-ploidy', '--ref-overlap', '--sample', sample]
    if regions_bed is not None:
        rtg_cmd += ['--bed-regions', regions_bed]
    call(rtg_cmd)

def classify_calls(rtg_bin, rtg_ref, somatic_vcf, truth_vcf, tumour_sample, temp_dir, regions_bed=None):
    eval_dir = join(temp_dir, "rtg-vcfeval")
    rtg_vcf_name = somatic_vcf.replace('.vcf', '.rtg.vcf')
    prepare_octopus_vcf_for_rtg(somatic_vcf, tumour_sample, rtg_vcf_name)
    run_rtg_vcfeval(rtg_bin, rtg_ref, rtg_vcf_name, truth_vcf, tumour_sample, eval_dir, regions_bed)
    remove_vcf(rtg_vcf_name)
    tp_vcf = join(temp_dir, "tp.vcf.gz")
    shutil.move(join(eval_dir, "tp.vcf.gz"), tp_vcf)
    fp_vcf = join(temp_dir, "fp.vcf.gz")
    shutil.move(join(eval_dir, "fp.vcf.gz"), fp_vcf)
    shutil.rmtree(eval_dir)
    return tp_vcf, fp_vcf

def subset(vcf_in_path, vcf_out_path, bed_regions):
    call(['bcftools', 'view', '-R', bed_regions, '-O', 'z', '-o', vcf_out_path, vcf_in_path])

def is_missing(x):
    if x == '.':
        return True
    x = float(x)
    return np.isnan(x)

def to_str(x, missing_value):
    if is_missing(x):
        return str(missing_value)
    else:
        return str(x)

def get_data(rec, features, n_samples, missing_value):
    result = [[] for _ in range(n_samples)]
    for feature in features:
        value = get_value(rec, feature)
        if type(value) == tuple:
            assert len(value) == n_samples
            result = [curr + [to_str(v, missing_value)] for curr, v in zip(result, value)]
        else:
            value_str = to_str(value, missing_value)
            for d in result:
                d.append(value_str)
    return result

def make_ranger_data(octopus_vcf_path, measures, is_tp, out, missing_value):
    vcf = VariantFile(octopus_vcf_path)
    n_samples = len(vcf.header.samples)
    n_records = 0
    with open(out, 'a') as ranger_dat:
        datwriter = csv.writer(ranger_dat, delimiter=' ')
        for rec in vcf:
            for row in get_data(rec, measures, n_samples, missing_value):
                if is_tp:
                    row.append('1')
                else:
                    row.append('0')
                datwriter.writerow(row)
                n_records += 1
    return n_records

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

def run_ranger_training(ranger, data_path, n_trees, min_node_size, threads, out, mtry=None):
    ranger_cmd = [ranger, '--file', data_path, '--depvarname', 'TP', '--probability',
                  '--ntree', str(n_trees), '--targetpartitionsize', str(min_node_size),
                  '--nthreads', str(threads), '--outprefix', out, '--write', '--verbose']
    if mtry is not None:
        ranger_cmd += ['--mtry', str(mtry)]
    call(ranger_cmd)

def main(options):
    tumour_samples = get_sample_names(options.truth)
    if len(tumour_samples) != 1:
        print("Truth VCF must contain one sample")
    tumour_sample = tumour_samples[0]
    if not exists(options.out):
        makedirs(options.out)
    master_data_path = join(options.out, options.prefix + ".dat")
    num_tps, num_fps = 0, 0
    for callset in options.variants:
        somatic_vcf_path = join(options.out, basename(callset).replace('.vcf', 'SOMATIC.tmp.vcf'))
        filter_somatic(callset, somatic_vcf_path)
        call(['tabix', somatic_vcf_path])
        tp_vcf_path, fp_vcf_path = classify_calls(options.rtg, options.sdf, somatic_vcf_path, options.truth, tumour_sample, options.out, options.regions)
        remove(somatic_vcf_path)
        remove(somatic_vcf_path + ".tbi")
        num_tps += make_ranger_data(tp_vcf_path, options.measures, True, master_data_path, options.missing_value)
        num_fps += make_ranger_data(fp_vcf_path, options.measures, False, master_data_path, options.missing_value)
        remove(tp_vcf_path)
        remove(fp_vcf_path)
    print("Number of TP examples: " + str(num_tps))
    print("Number of FP examples: " + str(num_fps))
    shuffle(master_data_path)
    ranger_header = ' '.join(options.measures + ['TP'])
    add_header(master_data_path, ranger_header)
    ranger_out_prefix = join(options.out, options.prefix)
    run_ranger_training(options.ranger, master_data_path, options.trees, options.min_node_size, options.threads, ranger_out_prefix, mtry=options.mtry)
    remove(ranger_out_prefix + ".confusion")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--variants',
                        type=str,
                        nargs='+',
                        required=True,
                        help='Octopus cancer calling model calls with CSR annotations')
    parser.add_argument('-T', '--regions',
                        type=str,
                        required=True,
                        help='BED files containing regions to use')
    parser.add_argument('--truth',
                        type=str,
                        required=True,
                        help='Truth VCF file')
    parser.add_argument('--measures',
                        type=str,
                        nargs='+',
                        required=True,
                        help='Measures to use for training')
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
                        required=True,
                        help='Output directory')
    parser.add_argument('--trees',
                        type=int,
                        default=300,
                        help='Number of trees to use in the random forest')
    parser.add_argument('--min_node_size',
                        type=int,
                        default=20,
                        help='Node size to stop growing trees, implicitly limiting tree depth')
    parser.add_argument('--mtry',
                        type=int,
                        default=None,
                        help='Number of variables to possibly split at in each node')
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
    main(parsed)
