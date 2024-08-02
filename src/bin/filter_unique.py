#!/usr/bin/env python

import sys
import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mixed_bam', help='reads mapped to mixed genome bam file')
parser.add_argument('-o', '--mixed_bam_unique', help='reads mapped to mixed genome, keep unique reads on host')

args = parser.parse_args()

mixed_bam = args.mixed_bam
mixed_bam_unique = args.mixed_bam_unique

chromosomes = (['chr' + str(i) for i in range(1, 20)] +
    ['chrX', 'chrY', 'chrM'] +  
    ['aav-is-standard-' + str(i) + '_woITR' for i in range(1, 7)])

with pysam.AlignmentFile(mixed_bam, "rb") as input_bam, \
    pysam.AlignmentFile(mixed_bam_unique, "wb", header=input_bam.header) as output_bam:
    for read in input_bam:
        read_name = read.query_name
        flag = read.flag
        mapq = read.mapping_quality
        chrom = read.reference_name
        ref_start = read.reference_start + 1 # convert to 1-based
        ref_end = read.reference_end
        if chrom in chromosomes or chrom == "chrV":
            if not (chrom in chromosomes and read.has_tag("XA")): # reads not mapped to multi locations of host
                output_bam.write(read)
