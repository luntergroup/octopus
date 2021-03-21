#!/usr/bin/env python3

import argparse
import pysam as ps
from pathlib import Path

def is_assigned_read(read):
    tags = dict(read.tags)
    return "HP" in tags and len(tags["HP"]) == 1

def filter_assigned_reads(in_bam_filename, out_bam_filename, region=None):
    in_bam = ps.AlignmentFile(in_bam_filename)
    out_bam = ps.AlignmentFile(out_bam_filename, 'wb', template=in_bam)
    contig, start, end = None, None, None
    if region is not None:
        if ':' in region:
            contig, rest = region.split(':')
            start, end = rest.split('-')
            start, end = start.replace(',', ''), end.replace(',', '')
            start, end = int(start), int(end)
        else:
            contig = region
    for read in in_bam.fetch(contig, start, end):
                if is_assigned_read(read):
                    out_bam.write(read)
    out_bam.close()
    ps.index(str(out_bam_filename))
            
def main(args):
    if args.in_bam == args.out_bam:
        print("--out-bam cannot be the same as --in-bam")
    else:
        filter_assigned_reads(args.in_bam, args.out_bam, args.region)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-I','--in-bam',
                        type=Path,
                        required=True,
                        help='Input Octopus realigned BAM')
    parser.add_argument('-O','--out-bam',
                        type=Path,
                        required=True,
                        help='Output BAM')
    parser.add_argument('-T','--region',
                        type=str,
                        required=False,
                        help='Region to filter')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
