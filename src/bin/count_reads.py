import sys
import pysam
import argparse
from collections import defaultdict
from cigar import Cigar
# from parse_cigar import get_read_clipped_positions, find_minimum_distance

infile, read_counts = sys.argv[1:]

reads_vector = set() # vector containing reads
reads_host = set() # host containing reads
reads_chimers = set() # any kinds of split reads, no filter
reads_vector_vector = set()
reads_host_host = set()
reads_vector_host = set()

exclude_flags = 0x4 | 0x100 | 0x800 # not primary reads
suppl_flag = 0x800
chromosomes = ['chr' + str(i) for i in range(1, 20)] + ['chrX', 'chrY', 'chrM']


def get_align_length(cigar_str):
    align_length = 0
    current_length = ""

    for char in cigar_str:
        if char.isdigit():
            current_length += char
        else:
            if char in ['M']: # D does not count
                align_length += int(current_length)
            current_length = ""

    return align_length

tot_pri_reads = 0
bamfile = pysam.AlignmentFile(infile, "rb")
for read in bamfile:
	read_name = read.query_name
	flag = read.flag
	chrom = read.reference_name
	strand = "-" if read.is_reverse else "+"
	query_len = read.query_length
	query_align_len = read.query_alignment_length
	cigar = read.cigarstring
	ref_start = read.reference_start + 1 # convert to 1-based

	if not flag & exclude_flags:
		if chrom in chromosomes or chrom == 'chrV':
			tot_pri_reads += 1
		if flag & 0x40: # first read
			read_name_se = read_name + "_1"
		if flag & 0x80: # second read
			read_name_se = read_name + "_2"
		
		if chrom == 'chrV':
			reads_vector.add(read_name_se)
		elif chrom in chromosomes:
			reads_host.add(read_name_se)

		if read.has_tag("SA"):
			sa = read.get_tag("SA")
			if len(sa.split(";")) == 2: # only 1 supplementary alignment
				# reads_chimers.add(read_name_se) # this include non canonical chromosome reads, e.g., chrUn_GL456368v1
				sa_primary_list = sa.split(";")[0].split(',') # only consider the primary/first alignment
				suppl_chrom, suppl_ref_start, suppl_strand, suppl_cigar = sa_primary_list[0:4]
				suppl_align_len = get_align_length(suppl_cigar)

				if chrom == 'chrV':
					if suppl_chrom == 'chrV':
						reads_chimers.add(read_name_se)
						reads_vector_vector.add(read_name_se)
						# print(f"{read_name_se},{flag},{chrom},{ref_start},{query_len},{query_align_len},{suppl_align_len},{cigar}; {sa} {sa_primary_list[0:4]}")
					elif suppl_chrom in chromosomes:
						reads_chimers.add(read_name_se)
						reads_host.add(read_name_se)
						reads_vector_host.add(read_name_se)
						# print(f"{read_name_se},{flag},{chrom},{ref_start},{query_len},{query_align_len},{suppl_align_len},{cigar}; {sa} {sa_primary_list[0:4]}")

				elif chrom in chromosomes:
					if suppl_chrom == 'chrV':
						reads_chimers.add(read_name_se)
						reads_vector.add(read_name_se)
						reads_vector_host.add(read_name_se)
						# print(f"{read_name_se},{flag},{chrom},{ref_start},{query_len},{query_align_len},{suppl_align_len},{cigar}; {sa} {sa_primary_list[0:4]}")
					elif suppl_chrom in chromosomes:
						reads_chimers.add(read_name_se)
						reads_host_host.add(read_name_se)
						# print(f"{read_name_se},{flag},{chrom},{ref_start},{query_len},{query_align_len},{suppl_align_len},{cigar}; {sa} {sa_primary_list[0:4]}")
			

num_reads_vector = len(reads_vector)
num_reads_host = len(reads_host)
num_reads_chimers = len(reads_chimers)
num_reads_vector_vector = len(reads_vector_vector)
num_reads_host_host = len(reads_host_host)
num_reads_vector_host = len(reads_vector_host)

with open(read_counts, 'w') as output_file:
	output_file.write(f"{tot_pri_reads},{num_reads_vector},{num_reads_host},{num_reads_chimers},{num_reads_vector_vector},{num_reads_host_host},{num_reads_vector_host}")
