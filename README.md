
# AAV Integration Project

## Project Overview

This project focuses on analyzing integration sites (IS) for Adeno-Associated Virus (AAV) using various bioinformatics tools and workflows. The project employs Snakemake workflows to handle different stages of data processing, including mapping, annotation, and coverage calculation.

## Directory Structure

```
**v1.0.0**
├── src
│   ├── Snakefile_annofinal.smk
│   ├── Snakefile_map.smk
│   ├── Snakefile_IS_spikeIn.smk
│   ├── Snakefile_cov.smk
│   ├── Snakefile_anno.smk
│   ├── Snakefile_IS.smk
│   ├── Snakefile_site_cov.smk
│   ├── run.sh
│   ├── bin
│   │   ├── plot_spikeInIS.html
│   │   ├── plotIS_human.rmd
│   │   ├── rm_collisions.py
│   │   ├── filter_and_label_unique.py
│   │   ├── IS_anno_v1.py
│   │   ├── deploy_index.R
│   │   ├── IS_anno_v1_human.py
│   │   ├── plotIS.rmd
│   │   ├── parse_cigar_v0.py
│   │   ├── plot_readstats.html
│   │   ├── circos_plot.rmd
│   │   ├── extract_readid.py
│   │   ├── plot_ISreads.rmd
│   │   ├── mk_md_report.py
│   │   ├── plot_readstats.rmd
│   │   ├── mk_index.py
│   │   ├── collapseIS.py
│   │   ├── plot_cov.R
│   │   ├── count_reads.py
│   │   ├── plot_chemericReads.rmd
│   │   ├── plotIS.html
│   │   ├── plot_spikeInIS.rmd
│   │   ├── extract_ids.py
│   │   ├── deploy_dash_cov_app_ms.py
│   │   ├── format_IS_w_sonication_reads.py
│   │   ├── sonication_sites.py
│   │   ├── parse_cigar_v1.py
│   │   ├── filter_host_proper_paired_reads.py
│   │   ├── filter_unique.py
│   │   ├── plot_trunc.rmd
│   │   ├── sonication_sites_v1.py
│   │   ├── parse_cigar.py
│   │   ├── IS_anno_final.py
│   │   ├── app_template_samples.py
│   │   └── IS_anno_final_human.py
│   ├── Snakefile_cleanIS.smk
│   ├── Snakefile.py
│   ├── Snakefile_plot.smk
│   ├── Snakefile_samples.py
│   ├── Snakefile_count.smk
│   ├── rmd
│   │   └── site
│   ├── Snakefile_ref.smk
│   ├── Snakefile_const.py
│   └── Snakefile_site.smk
├── data
│   ├── refs
│   ├── processed
│   ├── interim
│   │   └── .keep
│   └── raw
│       └── .keep
├── configs
│   ├── Cosmic_CancerGeneCensus_v98_GRCh38.tsv
│   ├── master_mapping.tsv
│   ├── AAV_integration_data_analyses_overview_original.xlsx
│   ├── AAV_integration_data_analyses_overview.xlsx
│   ├── index.Rmd
│   ├── master_mapping_anno.tsv
│   ├── master_mapping_rerun.tsv
│   └── aav_integration.json
├── docs
├── README.md
└── reqs
    ├── dash-cov-requirements.txt
    ├── aav_integration.conda.env.yaml
    ├── rsconnect.conda.env.yaml
    ├── samblaster.conda.env.yaml
    └── cat.conda.env.yaml
```

## Setup Instructions

### Prerequisites

Ensure you have the following software installed:
- [Conda](https://docs.conda.io/en/latest/miniconda.html)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Singularity](https://sylabs.io/guides/3.0/user-guide/)

### Installation

1. **Clone the repository**:
   ```sh
   git clone git@github.com:Sparktx-Data-Science/aav_integration.git
   cd aav_integration/v1.0.0
   ```

2. **Setup Conda Environment**:
   ```sh
   conda env create -f reqs/aav_integration.yaml
   conda activate aav_integration
   ```

3. **Setup Configurations**: 
  Edit the configuration files `aav_integration.json` in the `configs` directory as needed.

4. **Prepare Data**: 
   Ensure raw data is located in the appropriate directories as specified in the configuration.
   Symbol link raw data folder at /data/raw/shortRead/`run_id`

5. **Compile _site.yml file for deployment**:  
Store the following text in the file ./src/rmd/site/aav_integration_v1_0_0/_site.yml 
  ```sh
  name: "AAV integration"
  navbar:
    title: "AAV integration analysis report"
    left:
      - text: "About"
        href: about.html
    output:
        html_document:
          theme: darkly
  ```

6. **Fill out the experimental design file** 
   ```
   ./configs/AAV_integration_data_analyses_overview.xlsx
   ```
   
## Usage

### Running the Workflow

Execute the Snakemake workflow using the following command:

```sh
cd aav_integration/v1.0.0/
rm rmd/site/aav_integration_v1_0_0/index.Rmd
rm ../data/interim/flags/index_up.aav_integration_v1_0_0
rm ../data/interim/flags/rsconnect
bash run.sh
```

### Key Snakemake Rules

- **Mapping**: `Snakefile_map.smk`
  - Handles filtering, mapping reads to spike-in, vector, and host references, and removing duplicates.

- **Annotation**: `Snakefile_anno.smk`
  - Annotates integration sites, identifies nearest transcription start sites (TSS), and overlaps with genomic features.

- **Coverage Calculation**: `Snakefile_cov.smk`
  - Calculates coverage across different regions and generates visualizations.

- **Integration Sites Analysis**: `Snakefile_IS.smk`
  - Parses split reads, identifies integration sites, and performs sonication site analysis.

- **Spike-in Analysis**: `Snakefile_IS_spikeIn.smk`
  - Similar to IS analysis but specifically for spike-in controls.

### Example Commands

Run a specific Snakemake rule by modifying the run.sh script:
```sh
#!/bin/bash
# set -x

# 1. conda activate illumina_shortRead
# 2. prepare master_mapping.tsv in the ../configs folder
# 3. symbol link raw data folder at /data/raw/shortRead/

## ref__target map__target cov__target count__target IS_spikeIn__target IS__target anno__target cleanIS__target annofinal__target plot__target site__target

snakemake -s Snakefile.py \
        --use-conda --conda-frontend mamba --use-singularity \
        --singularity-args "-B /mnt/bfx-ops/:/mnt/bfx-ops/ -B /mnt/data/:/mnt/data/" \
        -j 60 \
        --configfile=../configs/aav_integration.json \
        --rerun-incomplete -p  \
        site__target 
```

### Generating Reports

The project includes several R Markdown (`.rmd`) scripts for generating detailed reports on integration sites, coverage, and other analyses. The reports are hosted on 
[Spark's Rstudio Connect server](https://connect.sparkds.io/aav_integration_v1_0_0/)

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on GitHub or contact chao.di@sparktx.com.
