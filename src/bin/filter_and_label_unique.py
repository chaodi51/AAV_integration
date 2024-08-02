#!/usr/bin/env python

import sys
import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mixed_bam', help='reads mapped to mixed genome bam file')
parser.add_argument('-o', '--mixed_bam_unique', help='reads mapped to mixed genome, keep unique reads on host')
parser.add_argument('-u', '--unique_label', help='read is labeled with unique or multi-mapping')

args = parser.parse_args()

mixed_bam = args.mixed_bam
mixed_bam_unique = args.mixed_bam_unique
unique_label = args.unique_label

# exclude_flags = 0x4 | 0x100 | 0x800 # not primary reads
chromosomes = (['chr' + str(i) for i in range(1, 20)] +
    ['chrX', 'chrY', 'chrM'] +  
    ['aav-is-standard-' + str(i) + '_woITR' for i in range(1, 7)])

with pysam.AlignmentFile(mixed_bam, "rb") as input_bam, \
    pysam.AlignmentFile(mixed_bam_unique, "wb", header=input_bam.header) as output_bam, \
    open(unique_label, 'w') as output_unique_label:
    for read in input_bam:
        read_name = read.query_name
        flag = read.flag
        read_name = read_name + '_1' if not flag & 0x80 else read_name + '_2'
        mapq = read.mapping_quality
        chrom = read.reference_name
        ref_start = read.reference_start + 1 # convert to 1-based
        ref_end = read.reference_end
        strand = "-" if read.is_reverse else "+"
        cigar = read.cigarstring

        if (chrom in chromosomes or chrom == "chrV"): # not only for primary alignments
            if chrom in chromosomes and read.has_tag("XA"): # reads mapped to multi locations of host
                output_unique_label.write(f"{read_name},{flag},{chrom},{ref_start},{strand},{cigar},host_multi\n")
            else: 
                output_bam.write(read) # reads map to chrV  or reads map to unique location of host
                if chrom in chromosomes:
                    output_unique_label.write(f"{read_name},{flag},{chrom},{ref_start},{strand},{cigar},host_uniq\n")
