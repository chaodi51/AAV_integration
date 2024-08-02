#!/usr/bin/env python

"""
In the v1 version, add uniqueness_lab for the IS reads 
"""

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
parser.add_argument('-u', '--uniq_label', help='host reads uniqueness labels')
parser.add_argument('-o', '--IS_tab', help='output table with IS')
parser.add_argument('-c', '--chimera_tab', help='output table with chimeras')

args = parser.parse_args()

file = args.input
uniq_label = args.uniq_label
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
## if primary alignment as the major alignment, it alway uses S clip, its supplementary alignment also use S clip
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
# print(get_ref_clip_positions("20S160M57S", 841))

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

def get_IS_position(read_name, clip_positions_on_read, host_align_clip_pos, host_strand):

    IS_pos = host_align_clip_pos # precise breakpoint or with Ins

    if read_name.endswith('_1'): # gap = read_clip_H - read_clip_V, read starts with Vector
        segment_gap = clip_positions_on_read[1] - clip_positions_on_read[0] 
        if host_strand == '+':
            if segment_gap < 0:
                IS_pos = host_align_clip_pos + abs(segment_gap)
        else:
            if segment_gap < 0:
                IS_pos = host_align_clip_pos - abs(segment_gap)
    elif read_name.endswith('_2'): # gap = read_clip_V - read_clip_H, read starts with Host
        segment_gap = clip_positions_on_read[0] - clip_positions_on_read[1]
        if host_strand == '+':
            if segment_gap < 0:
                IS_pos = host_align_clip_pos - abs(segment_gap)
        else:
            if segment_gap < 0:
                IS_pos = host_align_clip_pos + abs(segment_gap)

    return segment_gap, IS_pos

# Create a list to hold the results
IS_data = []
chimera_data = []
chromosomes = (['chr' + str(i) for i in range(1, 20)] +
    ['chrX', 'chrY', 'chrM'] +  
    ['aav-is-standard-' + str(i) + '_woITR' for i in range(1, 7)])

# retrive read uniqueness from earlier step
# "M05240:127:000000000-KVV5C:1:1101:10162:4721_2,2193,chr1,78585527,-,57M136H,host_uniq"
read_uniqueness = defaultdict(str)
with open(uniq_label, 'r') as uniq_label:
    for line in uniq_label:
        i = line.strip().split(',')
        raw_cigar = i[5]
        soft_clip_cigar = raw_cigar.replace("H", "S") if "H" in raw_cigar else raw_cigar
        read_uniqueness[i[0] + ',' + ','.join(i[2:5]) + ',' + soft_clip_cigar] = i[6]
        # print(','.join(i[0:5])+","+soft_clip_cigar+',',i[6])

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

    # get secondary alignments
    if read.has_tag("XA"):
        xa_tag = read.get_tag("XA")

    if read.has_tag("SA"):
        sa_tag = read.get_tag("SA")
        # only allow two parts when get split reads, so SA tag only has one alignment
        sa_list = sa_tag.rstrip(";").split(',')
        suppl_chrom, suppl_ref_start, suppl_strand, suppl_cigar = sa_list[0:4]
        sa_short = ','.join(sa_list[0:4])
        # read_len_cigar = get_read_length(suppl_cigar) # same as read_length
        
        # if uniqueness_lab == 'host_multi': # which must have XA
        #     # print(f"{read_name},{chrom},{ref_start},{strand},{cigar};--- {sa_tag}--- {xa_tag}")
        #     alternative_hits = xa_tag.rstrip(";").split(';')
        #     for hit in alternative_hits:
        #         components = hit.split(',')
        #         alt_chrom = components[0]
        #         alt_strand = components[1][0]  # First character indicates strand
        #         alt_ref_start = int(components[1][1:])  # Remaining characters indicate reference start
        #         alt_cigar = components[2]

                # read_clip_positions_alt = get_read_clipped_positions(alt_cigar, alt_strand)
                # ref_clip_positions_alt = get_ref_clip_positions(alt_cigar, alt_ref_start)

                # print(f"{alt_chrom},{alt_strand},{alt_ref_start},{alt_cigar}, {read_clip_positions_alt},{ref_clip_positions_alt}")
            # print(f"{chrom},{ref_start},{strand},{cigar}; {sa_tag}; {xa_tag}")


        read_clip_positions_pri = get_read_clipped_positions(cigar, strand)
        read_clip_positions_suppl = get_read_clipped_positions(suppl_cigar, suppl_strand)

        ref_clip_positions_pri = get_ref_clip_positions(cigar, ref_start)
        ref_clip_positions_suppl = get_ref_clip_positions(suppl_cigar, int(suppl_ref_start))

        min_clips_dist, min_clip_pair = find_minimum_distance(read_clip_positions_pri, read_clip_positions_suppl)
        segment_gap = min_clip_pair[0] - min_clip_pair[1]

        if strand == '+':
            ref_clip_pos_pri_major = ref_clip_positions_pri[read_clip_positions_pri.index(min_clip_pair[0])]
        else:
            ref_clip_pos_pri_major = ref_clip_positions_pri[::-1][read_clip_positions_pri.index(min_clip_pair[0])]

        if suppl_strand == '+':
            ref_clip_pos_suppl_major = ref_clip_positions_suppl[read_clip_positions_suppl.index(min_clip_pair[1])]
        else:
            ref_clip_pos_suppl_major = ref_clip_positions_suppl[::-1][read_clip_positions_suppl.index(min_clip_pair[1])]

        # if len(read_clip_positions_pri) > 1 or len(read_clip_positions_suppl) > 1:
        #     print(f"{chrom},{ref_start},{query_length},{strand},{cigar}; {sa} {read_clip_positions_pri},{read_clip_positions_suppl}; {min_clip_pair},{segment_gap}; {ref_clip_positions_pri},{ref_clip_positions_suppl},{ref_clip_pos_pri_major},{ref_clip_pos_suppl_major};")

        # if -30 < segment_gap < 30:
            # put chrV in the front as it is first sequenced
        if chrom in chromosomes and suppl_chrom in ['chrV']: # all IS reads
            vector_chrom = suppl_chrom
            vector_strand = suppl_strand
            vector_cigar = suppl_cigar
            vector_ref_start = suppl_ref_start
            host_chrom = chrom
            host_strand = strand
            host_cigar = cigar
            host_ref_start = ref_start
            clip_positions_on_read = min_clip_pair[::-1]
            vector_align_clip_pos = ref_clip_pos_suppl_major
            host_align_clip_pos = ref_clip_pos_pri_major
            # print(f"{vector_chrom}, {host_chrom}")

        else: # (chrom == 'chrV' and suppl_chrom in chromosomes) or (both from same species) # both IS and chimeric reads
            vector_chrom = chrom
            vector_strand = strand
            vector_cigar = cigar
            vector_ref_start = ref_start
            host_chrom = suppl_chrom
            host_strand = suppl_strand
            host_cigar = suppl_cigar
            host_ref_start = suppl_ref_start
            clip_positions_on_read = min_clip_pair
            vector_align_clip_pos = ref_clip_pos_pri_major
            host_align_clip_pos = ref_clip_pos_suppl_major
            # if vector_chrom == 'chrV' and host_chrom in chromosomes:
            #     print(f"{vector_chrom}, {host_chrom}")

        # label multi-mapper on host
        uniqueness_key = f"{read_name},{host_chrom},{host_ref_start},{host_strand},{host_cigar}"
        if uniqueness_key in read_uniqueness:
            uniqueness_lab = read_uniqueness[uniqueness_key]
        else:
            uniqueness_lab = "not_aval" # map to chrV
        # print(f"{uniqueness_key}, {uniqueness_lab}")


        segment_gap, IS_pos = get_IS_position(read_name, clip_positions_on_read, host_align_clip_pos, host_strand)

        data_row_IS = [read_name, uniqueness_lab, vector_chrom, vector_align_clip_pos, vector_strand, vector_cigar,
                    host_chrom, host_ref_start, IS_pos, host_strand, host_cigar, ';'.join(map(str, clip_positions_on_read)), segment_gap, query_seq]

        data_row_chimera = [read_name, vector_chrom, vector_align_clip_pos, vector_strand, vector_cigar,
                    host_chrom, host_ref_start, IS_pos, host_strand, host_cigar, ';'.join(map(str, clip_positions_on_read)), segment_gap, query_seq]

        if vector_chrom in ['chrV'] and host_chrom in chromosomes: ## IS reads
            # print(f"{vector_chrom}, {host_chrom}")
            if -30 <= segment_gap <= 30:
                IS_data.append(data_row_IS)
        else:
            chimera_data.append(data_row_chimera)            


# Create DataFrames
IS_df = pd.DataFrame(IS_data, columns=["read_name", "uniqueness_lab", "vector_chrom", "vector_align_clip_pos", "vector_strand",
                                       "vector_cigar", "host_chrom", "read_start", "IS_pos", "host_strand", "host_cigar",
                                       "clip_positions_on_read", "segment_gap", "IS_read"])

chimera_df = pd.DataFrame(chimera_data, columns=["read_name", "chrom", "ref_clip_pos_pri_major", "strand",
                                                 "cigar", "suppl_chrom", "read_start", "ref_clip_pos_suppl_major",
                                                 "suppl_strand", "suppl_cigar", "min_clip_pair", "segment_gap", "chimera_read"])

# Write DataFrames to CSV files
IS_df.to_csv(IS_tab, index=False)
chimera_df.to_csv(chimera_tab, index=False)