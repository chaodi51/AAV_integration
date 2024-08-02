rule IS__split_reads:
    input:
        # DATA + "processed/mixed_bams_rmdup/{platform}/{run}/{sample}__{ref}__mixed__rmdup.sorted.bam"
        DATA + "processed/mixed_bams_paired/{platform}/{run}/{sample}__{ref}__mixed_paired.sorted.bam"
        # DATA + "processed/mixed_bams_unique/{platform}/{run}/{sample}__{ref}__mixed__unique.sorted.bam"
    output:
        sam = DATA + "processed/mixed_bams_split_reads/{platform}/{run}/{sample}__{ref}__mixed.splitreads.sam",
    params:
        tmp = DATA + "processed/tmp/{sample}__{ref}__mixed.sam",
    container:
        "docker://mgibio/alignment_helper-cwl:2.2.1"
    threads: 8
    shell:
        '''
            mkdir -p "../data/processed/tmp/"
            # touch {params.tmp} 
            samtools sort -n -@ {threads} {input} | samtools view -h | samblaster -e --ignoreUnmated  -s {output.sam} \
            --maxSplitCount 2 --minIndelSize 30 --maxUnmappedBases 30 --minNonOverlap 20 -o {params.tmp} 
            rm -f {params.tmp}
        '''

rule IS__split_reads_parts:
    input:
        sam = DATA + "processed/mixed_bams_split_reads/{platform}/{run}/{sample}__{ref}__mixed.splitreads.sam",
    output:
        pri_sam = DATA + "processed/mixed_bams_split_reads_pri/{platform}/{run}/{sample}__{ref}__mixed.splitreads_pri.sam",
        # suppl_sam = DATA + "processed/mixed_bams_split_reads_suppl/{platform}/{run}/{sample}__{ref}__mixed.splitreads_suppl.sam",
    container:
        "docker://samesense/bwa-samtools:v1"
    threads: 8
    shell:
        '''
            samtools view -h -F 2048 {input} > {output.pri_sam}
        '''

rule IS__cigar_parse:
    input:
        pri_sam = DATA + "processed/mixed_bams_split_reads_pri/{platform}/{run}/{sample}__{ref}__mixed.splitreads_pri.sam",
        uniq_label = DATA + "processed/mixed_bams_host_unique_label/{platform}/{run}/{sample}__{ref}__mixed__host_unique_label.csv"
    output:
        IS_tab = INTs + "IS_results/IS_reads/{platform}/{run}/{sample}__{ref}__IS.csv",
        chimera_tab = INTs + "chimera_results/chimera_reads/{platform}/{run}/{sample}__{ref}__chimera.csv",
    shell:
        '''
            python bin/parse_cigar_v1.py -i {input.pri_sam} -u {input.uniq_label} -o {output.IS_tab} -c {output.chimera_tab}
        '''

# format to IS events with sonication abundance
rule IS__format_IS:
    input:
        IS = INTs + "IS_results/IS_reads/{platform}/{run}/{sample}__{ref}__IS.csv",
        bam = DATA + "processed/mixed_bams_rmdup/{platform}/{run}/{sample}__{ref}__mixed__rmdup.sorted.bam"
    output:
        INTs + "IS_results/IS_events/{platform}/{run}/{sample}__{ref}__IS.csv",
    params:
        lib_type = lambda wildcards: RUN_LIB_TYPE[wildcards.run]
    shell:
        '''
#         if echo {input} | grep -q '230818_VH00163_52_AACWHTHM5\|230925_M07542_0048_000000000-L6GCD\|231114_VH00163_66_AAF55KMM5\|\
# 231215_VH00163_75_AAF2M5TM5\|240112_VH00163_77_AAF5V37M5\|230816_M71053_29_000000000-LBGW4\|240209_VH00163_85_AAF5WC2M5\|\
# 240214_VH00163_87_AAF5W7NM5\|230817_LH00401_55_223JYKLT4_TES\|240329_VH00163_95_AAFCNC7M5\|\
# 240405_M05240_0149_000000000-L9FLR\|240328_VH00163_94_AAFJHJGM5'; 
        if [[ {params.lib_type} == "TES" ]];
        then
            echo -e 'TES library, both R1 and R2 can start with vector or host!!!\n'
            python bin/sonication_sites_v1.py -i {input.IS} -b {input.bam} -o {output} -m TES

        else
            python bin/sonication_sites_v1.py -i {input.IS} -b {input.bam} -o {output} -m slimPCR
        fi
        '''

## plot break points/IS positions on the vector
def IS_read_files(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        key = f"{platform}/{run}"
        if key == wc.run:
            for ref in SAMPLE_REFS[each]:
                ls.append(INTs + f"IS_results/IS_reads/{platform}/{run}/{sample}__{ref}__IS.csv")
    return ls

rule IS__plot_ISreads:
    input:
        IS_read_files
    output:
        INTs + "plots/{run}__ISreads.html", # run is platform/run
    params:
        IS_read = INTs + "IS_results/IS_reads/{run}/"
    shell:
        '''
        mkdir -p {INTs}plots/IS__plot_ISreads/{wildcards.run}/
        cp bin/plot_ISreads.rmd {INTs}plots/IS__plot_ISreads/{wildcards.run}/plot_ISreads.rmd
        
        Rscript -e '
            setwd("{INTs}plots/IS__plot_ISreads/{wildcards.run}/");
            rmarkdown::render("{INTs}plots/IS__plot_ISreads/{wildcards.run}/plot_ISreads.rmd", \
            output_format = "html_document", params = list(files = "{params.IS_read}"), \
            output_file = "{output}")'

        rm -r {INTs}plots/IS__plot_ISreads/{wildcards.run}/
        '''

## plot vector-vector break points positions on the vector
def chimeric_read_files(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        key = f"{platform}/{run}"
        if key == wc.run:
            for ref in SAMPLE_REFS[each]:
                ls.append(INTs + f"chimera_results/chimera_reads/{platform}/{run}/{sample}__{ref}__chimera.csv")
    return ls

rule IS__plot_chimericReads:
    input:
        chimeric_read_files,
        vector_anno = DATA + 'refs/vector/anno/CAG_EGFP_ITR_corrected.bed',
    output:
        INTs + "plots/{run}__chemericReads.html", # run is platform/run
    params:
        chimeric_read = INTs + "chimera_results/chimera_reads/{run}/"
    shell:
        '''
        mkdir -p {INTs}plots/IS__plot_chimericReads/{wildcards.run}/
        cp bin/plot_chemericReads.rmd {INTs}plots/IS__plot_chimericReads/{wildcards.run}/plot_chemericReads.rmd
        
        Rscript -e '
            setwd("{INTs}plots/IS__plot_chimericReads/{wildcards.run}/");
            rmarkdown::render("{INTs}plots/IS__plot_chimericReads/{wildcards.run}/plot_chemericReads.rmd", \
            output_format = "html_document", params = list(files = "{params.chimeric_read}",  anno = "{input.vector_anno}"), \
            output_file = "{output}")'

        rm -r {INTs}plots/IS__plot_chimericReads/{wildcards.run}/
        '''

# rule IS__plot_IS:
#     input:
#         IS_read = IS_read_files,
#         vector_anno = DATA + 'refs/vector/anno/CAG_EGFP_ITR_corrected.bed'
#     output:
#         INTs + "IS_results/IS_vector_plot/{run}__IS_cov_vector.html", # run is platform/run
#     # params:
#     #     vector_anno = DATA + 'refs/vector/anno/CAG_EGFP_ITR_corrected.bed'
#     shell:
#         '''
#         Rscript -e 'rmarkdown::render("bin/plot_trunc.rmd", \
#             output_format = "html_document", params = list(files = "{input.IS_read}", anno = "{input.vector_anno}"), \
#             output_file = "{output}")'
#         '''


# for middle step test
def collect_outputs(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        for ref in SAMPLE_REFS[each]:
            ls.append(INTs + f"IS_results/IS_reads/{platform}/{run}/{sample}__{ref}__IS.csv")
            ls.append(INTs + f"IS_results/IS_events/{platform}/{run}/{sample}__{ref}__IS.csv")
            # ls.append(INTs + f"chimera_results/chimera_reads/{platform}/{run}/{sample}__{ref}__chimera.csv")
    assert len(ls) >= 1
    return ls

def IS_plots_output(wc):
    ls = []
    for each in RUNS:
        platform, runid = each.split("--")
        run = platform + '/' + runid
        # ls.append(INTs + f"IS_results/IS_vector_plot/{platform}/{run}__IS_cov_vector.html")
        ls.append(INTs + f"plots/{run}__ISreads.html")
        ls.append(INTs + f"plots/{run}__chemericReads.html")
    assert len(ls) >= 1
    return ls

rule IS__target:
    input:
        collect_outputs,
        IS_plots_output

