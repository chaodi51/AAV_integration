#!/usr/bin/env python

import sys
import pysam
import numpy as np
import pandas as pd
import argparse
from cigar import Cigar
from collections import defaultdict
from itertools import combinations

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input sam file')
parser.add_argument('-o', '--IS_tab', help='output table with IS')
parser.add_argument('-c', '--chimera_tab', help='output table with chimeras')

args = parser.parse_args()

file = args.input
IS_tab = args.IS_tab
chimera_tab = args.chimera_tab

def get_read_length(cigar_str):
    read_length = 0
    current_length = ""

    for char in cigar_str:
        if char.isdigit():
            current_length += char
        else:
            if char in ['M', 'I', 'S', 'H']: # D does not count
                read_length += int(current_length)
            current_length = ""

    return read_length

## notes:
## in primary alignment as the major alignment, it alway uses S clip, its supplementary alignment also use S clip
## if the supplementary alignment as the major alignment, it use H clip, its supplementary alignment (the primary alignment) still use S clip

## get position of the split point relative to read start
def get_read_clipped_positions(cigar, strand):

    read_clip_pos = 1 # set to the read left most position
    read_clip_positions = []
    # e.g., [(20, 'H'), (20, 'M'), (20, 'S')]
    cigar_tuples = list(Cigar(cigar).items()) if strand == "+" else list(Cigar(cigar).items())[::-1]

    # Iterate through the CIGAR tuples
    for index, my_tuple in enumerate(cigar_tuples):
        cigar_op = my_tuple[1]
        cigar_len = my_tuple[0]
        # Check for the clip operations
        if cigar_op in ['S', 'H']:
            if index == 0:
                read_clip_pos += cigar_len # if 5' S/H, position relative to left most
            read_clip_positions.append(read_clip_pos)

        # Update the clip position based on the CIGAR operation
        if cigar_op in ['D', 'N']: # Operations that consume reference positions
            read_clip_pos += 0
        if cigar_op in ['M', 'I']:  # Operations that consume query positions
            read_clip_pos += cigar_len

    return read_clip_positions

## get position of the split point w.r.t the reference
def get_ref_clip_positions(cigar, ref_start):

    ref_clip_positions = []
    # e.g., [(20, 'H'), (20, 'M'), (20, 'S')]
    cigar_tuples = list(Cigar(cigar).items())

    # Iterate through the CIGAR tuples
    for index, my_tuple in enumerate(cigar_tuples):
        cigar_op = my_tuple[1]
        cigar_len = my_tuple[0]
        # Check for the clip operations
        if cigar_op in ['S', 'H']:
            if index == 0:
                ref_clip_pos = ref_start
            ref_clip_positions.append(ref_clip_pos)

        # Update the clip position based on the CIGAR operation
        if cigar_op in ['D', 'N']:
            ref_clip_pos += cigar_len
        if cigar_op in ['M', 'I']:
            if cigar_op == 'M':
                if index == 0:
                    ref_clip_pos = ref_start + cigar_len
                else:
                    ref_clip_pos += cigar_len
            if cigar_op == 'I':
                ref_clip_pos += 0

    return ref_clip_positions

# print(get_ref_clip_positions("115M97S", 71141490))

def find_minimum_distance(list1, list2):
    if not list1 or not list2:
        raise ValueError("Both lists must contain at least one element.")

    min_difference = float('inf')
    min_difference_pair = None

    for num1 in list1:
        for num2 in list2:
            difference = abs(num1 - num2)
            if difference < min_difference:
                min_difference = difference
                min_difference_pair = (num1, num2)

    return min_difference, min_difference_pair

# Create a list to hold the results
results_IS = []
results_chimera = []

# open the sam file with primary alignments
samfile = pysam.AlignmentFile(file, "r")
for read in samfile:
    read_name = read.query_name
    flag = read.flag
    chrom = read.reference_name
    ref_start = read.reference_start + 1 # convert to 1-based
    ref_end = read.reference_end
    query_start = read.query_alignment_start + 1  # convert to 1-based
    query_end = read.query_alignment_end

    query_length = read.query_length # included S not H
    query_seq = read.query_sequence
    strand = "-" if read.is_reverse else "+"
    mapq = read.mapping_quality
    cigar = read.cigarstring

    if read.has_tag("SA"):
        sa = read.get_tag("SA")
        sa_list = sa.rstrip("\;").split(',')
        suppl_chrom, suppl_ref_start, suppl_strand, suppl_cigar = sa_list[0:4]
        sa_short = ','.join(sa_list[0:4])
        # read_len_cigar = get_read_length(suppl_cigar) # same as read_length
        
        # print(f"{suppl_chrom},{suppl_ref_start},{suppl_strand},{suppl_cigar};")
        read_clip_positions_pri = get_read_clipped_positions(cigar, strand)
        read_clip_positions_suppl = get_read_clipped_positions(suppl_cigar, suppl_strand)

        ref_clip_positions_pri = get_ref_clip_positions(cigar, ref_start)
        ref_clip_positions_suppl = get_ref_clip_positions(suppl_cigar, int(suppl_ref_start))

        # if len(read_clip_positions_pri) > 1 or len(read_clip_positions_suppl) > 1:
        min_clips_dist, min_clip_pair = find_minimum_distance(read_clip_positions_pri, read_clip_positions_suppl)
        segment_gap = min_clip_pair[0] - min_clip_pair[1]
        ref_clip_pos_pri_major = ref_clip_positions_pri[read_clip_positions_pri.index(min_clip_pair[0])]
        ref_clip_pos_suppl_major = ref_clip_positions_suppl[read_clip_positions_suppl.index(min_clip_pair[1])]

        #combine other info to the IS

        
        # print(f"{chrom},{ref_start},{ref_end},{read_length},{strand},{cigar}; {sa} {read_clip_positions_pri},{read_clip_positions_suppl}; {min_clip_pair},{min_clips_dist}; {ref_clip_positions_pri},{ref_clip_positions_suppl};")
        
        #print(f"{chrom},{ref_start},{strand},{cigar},{query_length}; {sa_short}; {min_clip_pair},{segment_gap}; {ref_clip_pos_pri_major}; {ref_clip_pos_suppl_major};", file = output_file)
        
        # Append the results as a dictionary to the list
        results_dict = {
            "primary_align_read_id": read_name,
            "primary_align_chrom": chrom,
            "primary_align_start": ref_start,
            "primary_align_strand": strand,
            "primary_align_cigar": cigar,
            "query_length": query_length,
            "supplementary_alignment": sa_short,
            "clip_positions_on_read": min_clip_pair,
            "gap_between_segments": segment_gap,
            "primary_align_clip_pos": ref_clip_pos_pri_major,
            "supplementary_align_clip_pos": ref_clip_pos_suppl_major,
            "query_seq": query_seq
        }


        chromosomes = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']
        if (chrom in chromosomes and suppl_chrom == 'chrV') or (chrom == 'chrV' and suppl_chrom in chromosomes):
            # limit gap between [-30, 30]
            if -30 < segment_gap < 30:
                results_IS.append(results_dict)
        else:
            results_chimera.append(results_dict)

    else:
        pass

# Convert the list of dictionaries to a DataFrame
df_IS = pd.DataFrame(results_IS)
df_IS.to_csv(IS_tab, index=False) 

df_chimera = pd.DataFrame(results_chimera)
df_chimera.to_csv(chimera_tab, index=False) 