#!/bin/bash
# set -x

# 1. conda activate illumina_shortRead
# 2. prepare master_mapping.tsv in the ../configs folder
# 3. symbol link raw data folder at /data/raw/shortRead/
# 4. cp and edit _site.yml from other projects to SITE folder

## ref__target map__target cov__target count__target IS_spikeIn__target IS__target anno__target cleanIS__target annofinal__target plot__target site__target

snakemake -s Snakefile.py \
        --use-conda --conda-frontend mamba --use-singularity \
        --singularity-args "-B /mnt/bfx-ops/:/mnt/bfx-ops/ -B /mnt/data/:/mnt/data/" \
        -j 60 \
        --configfile=../configs/aav_integration.json \
        --rerun-incomplete -p -R annofinal__closest_host_TSS \
        site__target 
