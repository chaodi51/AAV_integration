#!/usr/bin/env python

import sys
import pandas as pd
import pysam
import argparse
from collections import defaultdict
from cigar import Cigar

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--IS_reads', help='input IS reads file')
parser.add_argument('-b', '--bam_file', help='all reads mapped to mixed genome')
parser.add_argument('-o', '--IS_events', help='output IS events file')
parser.add_argument('-m', '--lib_method', help='library method, i.e., TES or slimPCR')

args = parser.parse_args()

IS_reads = args.IS_reads
bam_file = args.bam_file
IS_events = args.IS_events
lib_method = str(args.lib_method)
print(lib_method)

IS_read_dict = defaultdict(str)
IS_paired_read_dict = defaultdict(str)
IS_event_dict = defaultdict(list)
sonication_reads_slimPCR = defaultdict(str)
sonication_reads_TES = defaultdict(str)

def cigar_match_position(cigar_string):
    cigar_list = [op for op in list(Cigar(cigar_string).items())]
    if sum(1 for _, op in cigar_list if op == 'S') == 0:
        M_pos = "full_match"
    elif sum(1 for _, op in cigar_list if op == 'S') == 1:
        if cigar_list[0][1] == 'S':
            M_pos = 'cigar_right'
        if cigar_list[-1][1] == 'S':
            M_pos = 'cigar_left'
    elif sum(1 for _, op in cigar_list if op == 'S') == 2:
        if cigar_list[0][0] < cigar_list[-1][0]:
            M_pos = 'cigar_left'
        else:
            M_pos = 'cigar_right'
    return M_pos

with open(IS_reads, 'r') as IS_reads:
    next(IS_reads)
    for line in IS_reads:
        l = line.strip().split(',')
        read_name = l[0]
        # read = "R1" if read_name.endswith('_1') else "R2"

        # didn't consider the insertion direction here
        # if consider the direction, should add congruent or inversed in the event_id (fragments mapped to different strands of chrV and host)
        event_id = l[5] + ',' + l[7]
        event_info = ";".join(l[:-1]) # remove read sequence
        # event_info = ";".join(l) # include read sequence
        IS_event_dict[event_id].append(event_info)
        IS_read_dict[read_name] = ""

        # IS'read paired read
        paired_read_name = read_name.replace("_1", "_2") if read_name.endswith('_1') else read_name.replace("_2", "_1")
        IS_paired_read_dict[paired_read_name] = ""

# exclude_flags = 0x4 | 0x100 | 0x800 # not primary reads
chromosomes = (['chr' + str(i) for i in range(1, 20)] +
    ['chrX', 'chrY', 'chrM'] +  
    ['aav-is-standard-' + str(i) + '_woITR' for i in range(1, 7)])

## pull out reads with sonication site
bam_file = pysam.AlignmentFile(bam_file, "rb")
for read in bam_file:
    read_name = read.query_name
    flag = read.flag
    read_name = read_name + '_1' if not flag & 0x80 else read_name + '_2'
    strand = "-" if read.is_reverse else "+"
    chrom = read.reference_name
    # cigar = read.cigarstring
    
    # query_start = read.query_alignment_start + 1  # convert to 1-based
    # query_end = read.query_alignment_end

    ##### important notes ######
    # for slimPCR method, sonication reads must be R2(flag & 0x80), and map to host
    # for TES method, sonication reads can be either R1 or R2, and map to host
    # sometimes, IS and sonication can be both on R1/R2, but the position can only be abtained from where the paired read start
    if (chrom in chromosomes) and read_name in IS_paired_read_dict:

        # M_pos = cigar_match_position(cigar)

        if strand == "+": # and (M_pos in {'cigar_left', "full_match"}):
            ref_start = read.reference_start + 1 # convert to 1-based
            site = ref_start # left most position
        elif strand == "-": # and (M_pos in {'cigar_right', "full_match"}):
            ref_end = read.reference_end # 0-based but points to one past the last aligned residue
            site = ref_end # right most position
        else:
            site = 'not_valid'

        if site != 'not_valid':
            if lib_method == 'slimPCR' and (flag & 0x80):
                sonication_read_slimPCR = ";".join([read_name, chrom, str(site)])
                # sonication_read_slimPCR = ";".join([read_name, chrom, str(site), query_seq])
                sonication_reads_slimPCR[read_name] = sonication_read_slimPCR
            elif lib_method == 'TES':
                sonication_read_TES = ";".join([read_name, chrom, str(site)])
                # sonication_read_TES = ";".join([read_name, chrom, str(site), query_seq])
                sonication_reads_TES[read_name] = sonication_read_TES     

        # elif strand == "-" and (M_pos in {'cigar_right', "full_match"}):
        #     site = ref_end # right most position
        #     if lib_method == 'slimPCR' and (flag & 0x80):
        #         sonication_read_slimPCR = ";".join([read_name, chrom, str(site)])
        #         sonication_reads_slimPCR[read_name] = sonication_read_slimPCR
        #     elif lib_method == 'TES':
        #         sonication_read_TES = ";".join([read_name, chrom, str(site)])
        #         sonication_reads_TES[read_name] = sonication_read_TES

        # print(f"{read_name}\t{flag}\t{strand}\t{chrom}\t{ref_start},{ref_end},{site}\t{cigar};")
# print(sonication_reads.values())

data = []
IS_info = 'IS_reads[read;vector_chrom;vector_align_clip_pos;vector_strand;vector_cigar;'\
            'host_chrom;read_start;IS_pos;host_strand;host_cigar;clip_positions_on_read;segment_gap]'

## matching sonication reads with IS reads
for event_id, event_info_list in IS_event_dict.items():
    chrom, IS_position = event_id.split(',')
    nReads = len(event_info_list)
    IS_sonication_reads = set() # unique sonication reads
    IS_sonication_sites = set()

    for each in event_info_list:
        IS_read_name = each.split(';')[0]
        soni_read_name = IS_read_name.replace("_1", "_2") if IS_read_name.endswith('_1') else IS_read_name.replace("_2", "_1")

        if lib_method == 'slimPCR' and soni_read_name in sonication_reads_slimPCR:
            IS_sonication_reads.add(sonication_reads_slimPCR[soni_read_name])
        elif lib_method == 'TES' and soni_read_name in sonication_reads_TES:
            IS_sonication_reads.add(sonication_reads_TES[soni_read_name])

    # get unique number of soniciate sites, different reads can have the same site
    for soni_read in IS_sonication_reads:
        _, soni_chr, soni_pos = soni_read.split(';')
        soni_site = soni_chr + ';' + soni_pos
        IS_sonication_sites.add(soni_site)

    nChops = len(IS_sonication_sites)

    if nReads >=1 and nChops >= 1: # keep IS events with >=1 sonication sites
        data.append([chrom, IS_position, nReads, nChops, event_info_list, list(IS_sonication_reads)])

df = pd.DataFrame(data, columns = ["chrom", "IS_pos", "nReads", "unique_sonications", IS_info, "sonication_sites"])
df.to_csv(IS_events, index=False)
