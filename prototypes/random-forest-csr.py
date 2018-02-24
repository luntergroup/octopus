#!/usr/bin/env python3

from pysam import VariantFile
import numpy as np
from sklearn.externals import joblib
import argparse
import os
from subprocess import call

# This script implements a prototype Random Forest filter for octopus calls.
# The input VCF is assumed to be an Octopus VCF file, annotated with at leas the measures
# given in the 'feastures' list below. This can be done with the octopus `--csr-train` command line option.

features = ['QUAL', 'AF', 'SB', 'GQ', 'DP', 'QD', 'FRF', 'CRF', 'GC', 'URF', 'MQ', 'MQ0', 'MQD', 'GT']

def parse_pileup(pileup):
    return [int(f) for f in pileup.split(',')]

def parse_pileups(pileups):
    return [parse_pileup(p) for p in pileups.split(':')]

def parse_bin(bin):
    strand_count_str, mq_hist_str, pileups = bin.split('|')
    strand_counts = tuple([int(count) for count in strand_count_str.split(',')])
    mq_hist = [int(q) for q in mq_hist_str.split(',')]
    return strand_counts, mq_hist, parse_pileups(pileups)

def parse_realignments(ra):
    return [parse_bin(bin) for bin in ra.replace('}','')[1:].split("{")]

def is_homozygous(gt):
    return all(a == 1 for a in gt)

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
        if feature == "RA":
            return parse_realignments(",".join(list(rec.info[feature])))
        else:
            result = rec.info[feature]
            if result == '.':
                return np.nan
            else:
                try:
                    return float(rec.info[feature])
                except ValueError:
                    return result
    elif feature in rec.samples[0]:
        if feature == 'GT':
            return int(is_homozygous(rec.samples[0][feature]))
        else:
            return int(rec.samples[0][feature])
    else:
        return np.nan

def read_vcf_record(rec, features, filtering=None):
    if filtering:
        if filtering == "PASS":
            if "PASS" in rec.filter:
                return [get_value(rec, feature) for feature in features]
        elif "PASS" not in rec.filter:
            return [get_value(rec, feature) for feature in features]
    else:
        return [get_value(rec, feature) for feature in features]

def read_vcf_data(vcf_path, features, max_records=None, include_pred=None, filtering=None):
    print("Reading " + str(len(features)) + " features from " + vcf_path)
    bcf_in = VariantFile(vcf_path)
    result = []
    for rec in bcf_in.fetch():
        if max_records is not None:
            if max_records == 0:
                break
            else:
                max_records -= 1
        x = read_vcf_record(rec, features, filtering)
        if include_pred is None or include_pred(x):
            result.append(x)
    print("Read " + str(len(result)) + " records")
    return result

def to_phred(ln_prob):
    ln10Div10 = .230258509299404568401799145468436420760110148862877297603
    return ln_prob / -ln10Div10

def filter_vcf_batch(batch, out_vcf, classifier, drop_info):
    X_batch = np.array([read_vcf_record(r, features) for r in batch])
    X_batch[np.isnan(X_batch)] = 0
    predictions = classifier.predict(X_batch)
    ln_probs = classifier.predict_log_proba(X_batch)
    num_passed = 0
    for i, rec in enumerate(batch):
        fields = str(rec).split('\t')
        if predictions[i]:
            fields[6] = "PASS"
            num_passed += 1
        else:
            fields[6] = "CSR"
        score_str = "{:0.2f}".format(to_phred(ln_probs[i][0]))
        if drop_info or len(fields[7]) == 0:
            fields[7] = "CSR=" + score_str
        else:
            fields[7] += ";CSR=" + score_str
        out_vcf.write("\t".join(fields))
    return num_passed

def filter_vcf(in_path, out_path, classifier, drop_info=False, batch_size=200000):
    vcf_in = VariantFile(in_path)
    new_header = vcf_in.header
    new_header.add_line('##FILTER=<ID=CSR,Description="CSR filtered">')
    new_header.add_line('##INFO=<ID=CSR,Number=1,Type=Float,Description="CSR score">')
    vcf_out = open(out_path, 'w')
    vcf_out.write(str(new_header))
    num_records_processed = 0
    num_passed = 0
    tick_size = 200000
    batch = []
    for rec in vcf_in:
        if len(batch) < batch_size:
            batch.append(rec)
        else:
            num_passed += filter_vcf_batch(batch, vcf_out, classifier, drop_info)
            num_records_processed += len(batch)
            if tick_size is not None and num_records_processed % tick_size == 0:
                print("Filtered " + str(num_records_processed) + " records (" + str(num_passed) + " PASS)")
            batch = [rec]
    if len(batch) > 0:
        num_passed += filter_vcf_batch(batch, vcf_out, classifier, drop_info)
        num_records_processed += len(batch)
    print("Done filtering " + str(num_records_processed) + " records")
    print("#records PASS: " + str(num_passed))
    print("#records FILTER: " + str(num_records_processed - num_passed))

def bgzip(src_vcf_path):
    call(['bgzip', src_vcf])
    call(['tabix', src_vcf + ".gz"])

def main(options):
    rf_model = joblib.load(options.model)
    compress_output = False
    if options.out_vcf.endswith(".gz"):
        options.out_vcf = options.out_vcf[:-3]
        compress_output = True
    filter_vcf(options.in_vcf, options.out_vcf, rf_model, options.drop_info)
    if compress_output:
        bgzip(options.out_vcf)
        os.remove(options.out_vcf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_vcf', type=str,
                        help='Input VCF')
    parser.add_argument('out_vcf', type=str,
                        help='Output VCF')
    parser.add_argument('model', type=str,
                        help='The SCIKIT random forest model')
    parser.add_argument('--drop_info', default=False,
                        action='store_true')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
