#!/usr/bin/env python

import sys
import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mixed_bam', help='reads mapped to mixed genome bam file')
parser.add_argument('-o', '--mixed_bam_paired', help='reads mapped to mixed genome, keep properly paired reads on host')

args = parser.parse_args()

mixed_bam = args.mixed_bam
mixed_bam_paired = args.mixed_bam_paired

chromosomes = (['chr' + str(i) for i in range(1, 20)] +
    ['chrX', 'chrY', 'chrM'] +  
    ['aav-is-standard-' + str(i) + '_woITR' for i in range(1, 7)])

with pysam.AlignmentFile(mixed_bam, "rb") as input_bam, \
    pysam.AlignmentFile(mixed_bam_paired, "wb", header=input_bam.header) as output_bam:
    for read in input_bam:
        read_name = read.query_name
        flag = read.flag
        read_name = read_name + '_1' if not flag & 0x80 else read_name + '_2'
        mapq = read.mapping_quality
        chrom = read.reference_name
        chrom_mate = read.next_reference_name
        ref_start = read.reference_start + 1 # convert to 1-based
        ref_end = read.reference_end
        strand = "-" if read.is_reverse else "+"
        cigar = read.cigarstring
        # print(f"{chrom}, {chrom_mate}")
        if (chrom in chromosomes or chrom == 'chrV') and (chrom_mate in chromosomes or chrom_mate == 'chrV'):
            if flag & 0x2 and chrom == chrom_mate: # properly paired on the same chrom
                output_bam.write(read)
            else: 
                if (chrom in chromosomes and chrom_mate == 'chrV') or (chrom == 'chrV' and chrom_mate in chromosomes):
                    if flag & 0x2:
                        output_bam.write(read)
