#!/usr/bin/env python

import sys
import pandas as pd
import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--IS_reads', help='input IS reads file')
parser.add_argument('-b', '--bam_file', help='all reads mapped to mixed genome')
parser.add_argument('-o', '--IS_events', help='output IS events file')

args = parser.parse_args()

IS_reads = args.IS_reads
bam_file = args.bam_file
IS_events = args.IS_events

IS_read_dict = defaultdict(str)
events = defaultdict(list)
sonication_sites = defaultdict(str)
sonication_reads = defaultdict(str)

with open(IS_reads, 'r') as IS_reads:
    next(IS_reads)
    for line in IS_reads:
        l = line.strip().split(',')
        read_name = l[0]
        # read = "R1" if read_name.endswith('_1') else "R2"

        # if don't consider the insertion direction
        # if consider the direction, should add congruent or inversed in the event_id (fragments mapped to different strands of chrV and host)
        event_id = ','.join(l[5:7])
        event_info = ";".join(l)
        # event_info = ";".join([read] + l[1:])
        events[event_id].append(event_info)

        read_name_ori = read_name.rstrip("_[12]")
        IS_read_dict[read_name_ori] = ""


exclude_flags = 0x4 | 0x100 | 0x800 # not primary reads
chromosomes = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']

bam_file = pysam.AlignmentFile(bam_file, "rb")
for read in bam_file:
    read_name = read.query_name
    flag = read.flag
    strand = "-" if read.is_reverse else "+"
    chrom = read.reference_name
    ref_start = read.reference_start + 1 # convert to 1-based
    ref_end = read.reference_end
    query_seq = read.query_sequence
    # query_start = read.query_alignment_start + 1  # convert to 1-based
    # query_end = read.query_alignment_end

    # sonication site must be on R2, and map to host
    if (not flag & exclude_flags) and (flag & 0x80) and (chrom in chromosomes): 
        if read_name in IS_read_dict:    
            if strand == "+":
                site = ref_start # left most position
            elif strand == "-":
                site = ref_end # right most position
            sonication_read = ";".join([read_name, chrom, str(site), query_seq])
            sonication_reads[read_name] = sonication_read
            sonication_site = ";".join([chrom, str(site)])
            sonication_sites[read_name] = sonication_site

data = []
IS_info = 'IS_reads[read;vector_chrom;vector_align_clip_pos;vector_strand;vector_cigar;'\
            'host_chrom;IS_pos;host_strand;host_cigar;clip_positions_on_read;segment_gap;read_seq]'

for event_id, event_info_list in events.items():
    chrom, IS_position = event_id.split(',')
    nReads = len(event_info_list)
    event_sonications = set()
    reads_sonications = set()
    for each in event_info_list:
        IS_read = each.split(';')[0].rstrip("_[12]")
        if IS_read in sonication_sites:
            event_sonications.add(sonication_sites[IS_read])
        if IS_read in sonication_reads:
            reads_sonications.add(sonication_reads[IS_read])

    nChops = len(event_sonications)
    if nReads >=2 and nChops >= 1: # keep IS events with >=2 sonication sites
        data.append({'chrom': chrom, 'IS_pos': IS_position, 'nReads': nReads, 'unique_sonications': nChops, IS_info: event_info_list, "sonication_sites": list(reads_sonications)})

df = pd.DataFrame(data)
df.to_csv(IS_events, index=False)
