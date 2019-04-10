#!/usr/bin/env python3

import argparse
import pysam as ps
from pathlib import Path

def get_haplotype_id(read):
    try:
        return int(read.get_tag('hi'))
    except KeyError:
        return None

def main(options):
    in_bam = ps.AlignmentFile(options.bam)
    out_bam_prefix = Path(options.output)
    variant_out_bams = [ps.AlignmentFile(str(out_bam_prefix) + '_0.bam', 'wb', template=in_bam)]
    ref_out_bam = ps.AlignmentFile(str(out_bam_prefix) + '.bam', 'wb', template=in_bam)
    for read in in_bam:
        haplotype_id = get_haplotype_id(read)
        if haplotype_id is None:
            ref_out_bam.write(read)
        else:
            if len(variant_out_bams) < haplotype_id + 1:
                for n in range(len(variant_out_bams), haplotype_id + 1):
                    variant_out_bams.append(ps.AlignmentFile(str(out_bam_prefix) + '_' + str(n) + '.bam', 'wb', template=in_bam))
            variant_out_bams[haplotype_id].write(read)
    for bam in [ref_out_bam] + variant_out_bams:
        bam.close()
        ps.index(bam.filename)

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
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
