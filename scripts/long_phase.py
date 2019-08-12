#!/usr/bin/env python3

import argparse
import os
import os.path
import subprocess as sp

def vcf_index_exists(vcf_fname):
    return os.path.exists(vcf_fname + '.tbi')

def remove_vcf_index(vcf_fname):
    os.remove(vcf_fname + '.tbi')

def remove_vcf(vcf_fname):
    os.remove(vcf_fname)
    if vcf_index_exists(vcf_fname):
        remove_vcf_index(vcf_fname)

def main(options):
    octopus_cmd = [options.octopus,
                   '-R', options.reference,
                   '-I', options.reads]
    if options.regions is not None:
        octopus_cmd += ['--regions'] + options.regions
    if options.regions_file is not None:
        octopus_cmd += ['--regions-file', options.regions_file]
    if options.threads is not None:
        octopus_cmd += ['--threads', options.threads]
    if options.forest is not None:
        octopus_cmd += ['--forest', options.forest]
    if options.config is not None:
        octopus_cmd += ['--config', options.config]
    
    # first pass
    tmp_vcf = options.output.replace('.vcf', '.tmp.vcf')
    first_pass_cmd = octopus_cmd + ['-o', tmp_vcf]
    sp.call(first_pass_cmd)
    
    # second pass
    second_pass_cmd = octopus_cmd
    second_pass_cmd += ['--disable-denovo-variant-discover']
    second_pass_cmd += ['--source-candidates', tmp_vcf]
    second_pass_cmd += ['--lagging-level', 'AGGRESSIVE', '--backtrack-level', 'AGGRESSIVE']
    second_pass_cmd += ['-o', options.output]
    if options.bamout is not None:
        second_pass_cmd += ['--bamout', options.bamout]
        if options.bamout_type:
            second_pass_cmd += ['--bamout-type', options.bamout_type]
    sp.call(second_pass_cmd)

    # cleanup
    remove_vcf(tmp_vcf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('octopus',
                        type=str,
                        help='Octopus binary')
    parser.add_argument('-R','--reference',
                        type=str,
                        required=True,
                        help='Reference file in FASTA format')
    parser.add_argument('-I', '--reads',
                        type=str,
                        required=True,
                        help='Aligned reads in BAM format')
    parser.add_argument('-o', '--output',
                        type=str,
                        required=True,
                        help='Where to write the final output VCF')
    parser.add_argument('-T', '--regions',
                        nargs='+',
                        type=str,
                        required=False,
                        help='BED file containing regions to call/phase')
    parser.add_argument('-t', '--regions-file',
                        type=str,
                        required=False,
                        help='BED file containing regions to call/phase')
    parser.add_argument('--threads',
                        type=int,
                        required=False,
                        help='Number of threads to use for analysis')
    parser.add_argument('--forest',
                        type=str,
                        required=False,
                        help='Octopus random forest for filtering calls')
    parser.add_argument('--bamout',
                        type=str,
                        required=False,
                        help='Output realigned evidence BAM')
    parser.add_argument('--bamout-type',
                        type=str,
                        required=False,
                        help='Type of realigned evidence BAM to make')
    parser.add_argument('--config',
                        type=str,
                        required=False,
                        help='Config file to use for calling')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)