#!/usr/bin/env python3

import argparse
import os
import pysam as ps
from copy import copy

def copy_cigar(read, pos, length):
    result = []
    for op_id, op_length in read.cigar:
        if pos < op_length:
            op_length -= pos
            result.append((op_id, min(length, op_length)))
            if op_id not in [2, 5]: # deletion or hard clip 
                if length <= op_length:
                    break
                length -= op_length
                if pos > 0:
                    pos = 0
        elif op_id not in [2, 5]: # deletion or hard clip 
            pos -= op_length
    return result

def reference_len(cigar):
    return sum(op_length for op_id, op_length in cigar if op_id in [0, 2, 6, 7, 8])

def split_read(read, chunk_len):
    read_len = len(read.query_sequence)
    if chunk_len >= read_len:
        return [read]
    result = []
    chunk_reference_start = read.reference_start
    for idx in range(int(read_len / chunk_len)):
        sequence_offset = idx * chunk_len
        chunk_cigar = copy_cigar(read, sequence_offset, chunk_len)
        chunk_ref_len = reference_len(chunk_cigar)
        if chunk_ref_len == 0:
            continue
        chunk = copy(read)
        chunk.reference_start = chunk_reference_start
        chunk_reference_start += chunk_ref_len
        chunk.query_name += '_' + str(idx)
        chunk.tags += [('BX', read.query_name)]
        chunk.cigar = chunk_cigar
        chunk_seq = chunk.query_sequence[sequence_offset : sequence_offset + chunk_len]
        chunk_quals = chunk.query_qualities[sequence_offset : sequence_offset + chunk_len]
        chunk.query_sequence, chunk.query_qualities = chunk_seq, chunk_quals
        result.append(chunk)
    return result

def main(options):
    in_bam = ps.AlignmentFile(options.input)
    tmp_bam = ps.AlignmentFile(options.output + '.unsorted.bam', 'wb', template=in_bam)
    for read in in_bam:
        read.query_qualities = [options.base_quality for _ in range(len(read.query_sequence))]
        for chunk in split_read(read, options.chunk_length):
            tmp_bam.write(chunk)
    in_bam.close()
    tmp_bam.close()
    ps.sort("-o", options.output, tmp_bam.filename)
    os.remove(tmp_bam.filename)
    ps.index(options.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        type=str,
                        help='Input PacBio BAM file')
    parser.add_argument('output',
                        type=str,
                        help='Output BAM file')
    parser.add_argument('-l', '--chunk-length',
                        type=int,
                        default=300,
                        help='Split PacBio reads into chunks of this length')
    parser.add_argument('-q', '--base-quality',
                        type=int,
                        default=30,
                        help='Reset all base qualities to this')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
