#!/usr/bin/env python

import sys
import csv
import pandas as pd
import pysam
import argparse
import ast
from collections import defaultdict
# Increase the maximum field size limit
csv.field_size_limit(sys.maxsize)

events_file, host_tss, host_ingene, \
enh_ingene, cpg_ingene, rmsk_anno, output_file = sys.argv[1:]

events_dict = defaultdict(list)

with open(host_tss, 'r') as host_tss_file:
    for line in host_tss_file:
        l = line.strip().split('\t')
        event = l[0]+ ';' + l[2]
        host_tss = l[6] + "@" + str(l[9])
        events_dict[event] = [host_tss]

with open(host_ingene, 'r') as host_ingene_file:
    for line in host_ingene_file:
        l = line.strip().split('\t')
        event = l[0]+ ';' + l[2]
        host_ingene = l[3]
        events_dict[event].append(host_ingene)

with open(enh_ingene, 'r') as enh_ingene_file:
    for line in enh_ingene_file:
        l = line.strip().split('\t')
        event = l[0]+ ';' + l[2]
        enh_ingene = l[3]
        events_dict[event].append(enh_ingene)

with open(cpg_ingene, 'r') as cpg_ingene_file:
    for line in cpg_ingene_file:
        l = line.strip().split('\t')
        event = l[0]+ ';' + l[2]
        cpg_ingene = l[3]
        events_dict[event].append(cpg_ingene)

with open(rmsk_anno, 'r') as rmsk_anno_file:
    for line in rmsk_anno_file:
        l = line.strip().split('\t')
        event = l[0]+ ';' + l[2]
        rmsk_anno = l[3]
        events_dict[event].append(rmsk_anno)

sample, vector = events_file.split('/')[-1].rstrip('__IS.csv').split('__')
with open(events_file, 'r') as events_file, open(output_file, 'w', newline='') as outfile:
    csv_reader = csv.reader(events_file)
    header = next(csv_reader)
    # reorder columns
    new_header = ["sample", "vector", "event_lab", "IS_origin", "chrom", "IS_pos", "vector_align_clip_pos", "nReads", "unique_sonications", 
                  "host_closest_TSS_info", "host_closest_TSS_dist",
                  "host_inGene", "enh_ingene", "cpg_ingene", "rmsk(name/class/family)",
                  "IS_reads[read;vector_chrom;vector_align_clip_pos;vector_strand;vector_cigar;host_chrom;IS_pos;host_strand;host_cigar;clip_positions_on_read;segment_gap]",
                  "sonication_sites"]
    writer = csv.writer(outfile)
    writer.writerow(new_header)

    for l in csv_reader:
        chrom, IS_pos = l[2], l[3]
        id = f"{chrom};{IS_pos}"
        IS_reads = ast.literal_eval(l[6])
        
        vector_align_clip_pos = str([i.split(";")[3] for i in IS_reads])
        sonication_sites = l[7]

        if id in events_dict:
            host_tss_info, host_tss_dist  = events_dict[id][0].split("@")
            host_inGene, enh_ingene, cpg_ingene, rmsk_anno =  events_dict[id][1:5]
            
            row = [sample, vector] + l[0:4] + [vector_align_clip_pos] + l[4:6] + [host_tss_info, host_tss_dist, host_inGene, enh_ingene, cpg_ingene, rmsk_anno, str(IS_reads), sonication_sites]
            writer.writerow(row)

