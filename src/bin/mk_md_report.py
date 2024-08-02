"""Make rmd report for run"""
import sys
import csv
import os
import pandas as pd
from datetime import date

# these arguments have single file as input
(site_dir,          # SITE
run,                # run name
sample_ref,         # sample ref matchup file: MAP_FILE
merge_stats,        # merged reads stats from the pre-preprocess report
qc,                 # multiqc report from preprocess (.html)
map_stats,          # mapping stats (.tsv)
merged_IS,
merged_IS_final,
merged_IS_spikeIn,
readStats_plot,
ISreads_plot,
spikeInIS_plot,
IS_plot_per_sample,
IS_plot_combine_replicates
) = sys.argv[1:15]

out_rmd = sys.argv[-1]            # output run.rmd file

"""
these input arguments have multiple files each:
IS_anno,        # IS events with annotation for each sample (.csv)
IS_spikeIn,     # IS events identified from spike-Ins
cov_ip_file,    # coverage app ip files ([.cov.ip])
"""

merge_stats = pd.read_csv(merge_stats, sep = "\t")
merge_stats_tab = merge_stats.to_markdown(index = False, tablefmt = "pipe")

map_stats_csv = pd.read_csv(map_stats, sep = "\t")
map_stats_tab = map_stats_csv.to_markdown(index = False, tablefmt = "pipe")

sample_ref = pd.read_csv(sample_ref, sep = "\t")
sample_ref = (sample_ref[sample_ref["run"] == run].rename(columns = {"mapLsSemiDelim": "references"}))
sample_ref_tab = sample_ref.to_markdown(index = False, tablefmt = "pipe")

# mapping stats table
file_name = map_stats.split('/')[-1]
report_name = file_name.split('.')[0]
map_stats_tab_link = f"* [{report_name}]({file_name})\n"
os.system(f'cp {map_stats} {site_dir}/{file_name}')

# readStats_plot
file_name = readStats_plot.split('/')[-1]
report_name = file_name.split('.')[0]
readStats_plot_link = f"* [{report_name}]({file_name})\n"
os.system(f'cp {readStats_plot} {site_dir}/{file_name}')

# ISreads_plot
file_name = ISreads_plot.split('/')[-1]
report_name = file_name.split('.')[0]
ISreads_plot_link = f"* [{report_name}]({file_name})\n"
os.system(f'cp {ISreads_plot} {site_dir}/{file_name}')

# spikeInIS_plot
file_name = spikeInIS_plot.split('/')[-1]
report_name = file_name.split('.')[0]
spikeInIS_plot_link = f"* [{report_name}]({file_name})\n"
os.system(f'cp {spikeInIS_plot} {site_dir}/{file_name}')


# IS events with annotation (.csv)
IS_tables = [i for i in sys.argv if 'IS_events_anno' in i and 'spikeIn' not in i and 'final' not in i]
IS_anno_tab = ""
for file in IS_tables:
    run_id = file.split('/')[-2]
    file_name = file.split('/')[-1]
    report_name = file_name.split('.')[0]
    full_name = run_id + '__' + file_name
    
    IS_anno_tab += f"* [{report_name}]({full_name})\n"
    # os.system(f'cp {file} {site_dir}/{full_name}')

# merged IS table
file_name = merged_IS.split('/')[-1]
report_name = file_name.split('.')[0]
merged_IS_tab = f"* [{report_name}]({file_name})\n"
os.system(f'cp {merged_IS} {site_dir}/{file_name}')

# merged final IS table
file_name = merged_IS_final.split('/')[-1]
report_name = file_name.split('.')[0]
merged_IS_final_tab = f"* [{report_name}]({file_name})\n"
os.system(f'cp {merged_IS_final} {site_dir}/{file_name}')

# IS plot per sample
file_name = IS_plot_per_sample.split('/')[-1]
report_name = file_name.split('.')[0]
IS_plot_per_sample_link = f"* [{report_name}]({file_name})\n"
os.system(f'cp {IS_plot_per_sample} {site_dir}/{file_name}')

# IS plot combine replicates
file_name = IS_plot_combine_replicates.split('/')[-1]
report_name = file_name.split('.')[0]
IS_plot_combine_replicates_link = f"* [{report_name}]({file_name})\n"
os.system(f'cp {IS_plot_combine_replicates} {site_dir}/{file_name}')

# merged_IS_spikeIn
file_name = merged_IS_spikeIn.split('/')[-1]
report_name = file_name.split('.')[0]
merged_IS_spikeIn_tab = f"* [{report_name}]({file_name})\n"
os.system(f'cp {merged_IS_spikeIn} {site_dir}/{file_name}')

# IS spikeIn events (.csv)
IS_spikeIn_tables = [i for i in sys.argv if 'spikeIn_IS_events_anno' in i]
IS_spikeIn_tab = ""
for file in IS_spikeIn_tables:
    run_id = file.split('/')[-2]
    file_name = file.split('/')[-1]
    report_name = file_name.split('.')[0]
    full_name = run_id + '__' + report_name + '_spikeIn.csv'
    
    IS_spikeIn_tab += f"* [{report_name}]({full_name})\n"
    # os.system(f'cp {file} {site_dir}/{full_name}')

# coverage app ip files (.cov.ip)
cov_ip_file = [i for i in sys.argv if 'cov.ip' in i]
cov_app = ""
for ip_file in cov_ip_file:
    file_name = ip_file.split('/')[-1]
    report_name = file_name.split('--')[2] ## vector name
    with open(ip_file) as f:
        aav_cov_url = f.readline().strip()
        # this cannot be empty
        # if it is, app deployment have failed
        assert aav_cov_url
    cov_app += f"* [Interactive coverage of {report_name}]({aav_cov_url})\n"


template = f"""---
title: "{run}"
---

```{{r setup, include = FALSE}}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(ggplot2)
library(dplyr)
```

### Sequencing details
* Run name: {run}
* System: illumina mseq/nseq/novaseq
* Service: AAV integration
* Producer: Spark Data Science

### Sample reference matchup
{sample_ref_tab}

### Fastqc report
* [Multiqc report]({qc.split('/')[-1]})

### Merge reads stats
{merge_stats_tab}

### Read count table
{map_stats_tab}

Download table above:
{map_stats_tab_link}

### Read stats plots
{readStats_plot_link}

### IS reads composition plots
For a single chimeric read, the composition of AAV and host sequences can be:  
* gap == 0 ~ 'Precise Breakpoint', AAV directly connected to host genome.  
* gap > 0 ~ 'Insertion', inserted nucleotides between AAV and host genome.  
* gap < 0 ~ 'Micro-homology', overlapped nucleotides between AAV and host genome.  
Here we restrained the gap between (-30, 30)  

{ISreads_plot_link}

### IS events with annotations
Table with all IS events, read abundance, closest genes, CpG, enhancer, repeat features etc.

Merged IS table:
{merged_IS_tab}

### Merge IS within 10bp and remove collisions among unrelated samples
{merged_IS_final_tab}

### Feature plots for Merge IS events (per sample)
{IS_plot_per_sample_link}

### Feature plots for Merge IS events (combine replicates)
{IS_plot_combine_replicates_link}

### IS events detected in spikeIns

Merged IS table:
{merged_IS_spikeIn_tab}

### Plots for IS detected using spikeIns/standards

{spikeInIS_plot_link}

### Interactive coverage plot for each vector

{cov_app}

"""

with open(out_rmd, "w") as fout:
    print(template, file = fout)

