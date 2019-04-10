#!/usr/bin/env python3

import argparse
import pysam as ps
from pathlib import Path
from os.path import isfile
from os import rename, remove

def get_haplotype_id(read):
    try:
        return int(read.get_tag('hi'))
    except KeyError:
        return None

def bam_index_exists(bam_fname):
    return isfile(str(bam_fname.with_suffix(bam_fname.suffix + '.bai')))

def index_bam(bam_fname, sort_if_needed=False):
    ps.index(str(bam_fname))
    if not bam_index_exists(bam_fname):
        sorted_bam_fname = str(bam_fname).replace('.bam', '.sorted.tmp.bam')
        ps.sort('-o', sorted_bam_fname, str(bam_fname))
        rename(sorted_bam_fname, bam_fname)
        ps.index(str(bam_fname))

def main(options):
    if options.assigned_only and options.merge_reference:
        print('Options --assigned_only and --merge_reference are mutually exclusive')
        exit(0)
    if not isfile(options.bam):
        print('Input realigned BAM ' + options.bam + ' does not exist')
        exit(1)
    
    in_bam = ps.AlignmentFile(options.bam)
    out_bam_prefix = Path(options.output)
    variant_out_bams = [ps.AlignmentFile(str(out_bam_prefix) + '_0.bam', 'wb', template=in_bam)]
    ref_out_bam = None if options.assigned_only else ps.AlignmentFile(str(out_bam_prefix) + '.bam', 'wb', template=in_bam)
    
    for read in in_bam:
        haplotype_id = get_haplotype_id(read)
        if haplotype_id is None:
            if ref_out_bam is not None:
                ref_out_bam.write(read)
        else:
            if len(variant_out_bams) < haplotype_id + 1:
                for n in range(len(variant_out_bams), haplotype_id + 1):
                    variant_out_bams.append(ps.AlignmentFile(str(out_bam_prefix) + '_' + str(n) + '.bam', 'wb', template=in_bam))
            variant_out_bams[haplotype_id].write(read)
    
    if ref_out_bam is not None and options.merge_reference:
        ref_out_bam.close()
        merged_bam_fname = options.output + '.merged.tmp.bam'
        variant_out_bams[-1].close()
        unassigned_variant_bam_fname = variant_out_bams[-1].filename.decode("utf-8")
        ps.merge(merged_bam_fname, ref_out_bam.filename, unassigned_variant_bam_fname)
        rename(merged_bam_fname, unassigned_variant_bam_fname)
        remove(ref_out_bam.filename)
        ref_out_bam = None
        
    for bam in [ref_out_bam] + variant_out_bams:
        if bam is not None:
            bam.close()
            index_bam(Path(bam.filename.decode("utf-8")), sort_if_needed=options.sort)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bam',
                        type=str,
                        required=True,
                        help='Octopus realigned BAM file to split')
    parser.add_argument('-o','--output',
                        type=str,
                        required=True,
                        help='Output prefix')
    parser.add_argument('-A', '--assigned_only',
                        default=False,
                        action='store_true',
                        help='Do not output unassigned reads')
    parser.add_argument('-M', '--merge_reference',
                        default=False,
                        action='store_true',
                        help='Merge reference reads into unassigned variant reads')
    parser.add_argument('-S', '--sort',
                        default=False,
                        action='store_true',
                        help='Sort output by position if required')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
