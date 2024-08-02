include: "Snakefile_site_cov.smk"

rule site__init_rsconnect:
    output:
        touch(INTs + 'flags/rsconnect')
    conda:
        REQS + 'rsconnect.conda.env.yaml'
    shell:
        '''
        rsconnect add --api-key za4RPQdH58QuGHPTZbHvU087nF4TlkD3 \
            --server https://connect.sparkds.io/ --name rsconnect
        '''

# def mk_all_multiqc(wc):
#     ls = []
#     for pr in get_platform_runs():
#         html = QC + f'multiqc/{pr}.multiqc_report.html'
#         ls.append(html)
#     return ls

rule site__cp_site:
    input:
        html = QC + 'multiqc/shortRead--{run}.multiqc_report.html'
    output:
        SITE + 'shortRead--{run}.multiqc_report.html'
    shell:
        "cp {input} {output}"

def IS_anno_output(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        if run == wc.run:
            for ref in SAMPLE_REFS[each]:
                ls.append(INTs + f"IS_results/IS_events_anno/{platform}/{run}/{sample}__{ref}__IS.csv")
    assert len(ls) >= 1
    return ls

def IS_spikeIn_output(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        if run == wc.run:
            ls.append(INTs + f"spikeIn_IS_results/spikeIn_IS_events_anno/{platform}/{run}/{sample}__IS.csv")
    assert len(ls) >= 1
    return ls

rule site__run_rmd:
    input:
        merge_tab = MERGE_TAB + "{run}__merge_table.tsv",
        qc = SITE + 'shortRead--{run}.multiqc_report.html',
        map_stats = INTs + "shortRead/{run}__mapped_table.tsv",
        merged_IS = INTs + "IS_results/IS_events_comb/shortRead/{run}__IS.anno.csv",
        merged_IS_final = INTs + "IS_results_final/IS_events_anno/shortRead/{run}__IS_final.anno.csv",
        merged_IS_spikeIn = INTs + "spikeIn_IS_results/spikeIn_IS_events_comb/shortRead/{run}__spikeIn_IS.csv",
        readStats_plot = INTs + 'plots/shortRead/{run}__readStats.html',
        ISreads_plot = INTs + "plots/shortRead/{run}__ISreads.html",
        spikeInIS_plot = INTs + "plots/shortRead/{run}__spikeInIS.html",
        IS_plot_per_sample = INTs + 'plots/shortRead/{run}__IS_per_sample.html',
        IS_plot_combine_replicates = INTs + 'plots/shortRead/{run}__IS_combine_replicates.html',
        IS = IS_anno_output,
        IS_spikeIn = IS_spikeIn_output,
        cov = mk_cov_vis_ls_for_run
    output:
        SITE + 'shortRead--{run}.Rmd'
    params:
        sample_ref = MAP_FILE
    shell:
        '''
        python bin/mk_md_report.py {SITE} {wildcards.run} {params.sample_ref} \
            {input.merge_tab} {input.qc} {input.map_stats} {input.merged_IS} {input.merged_IS_final} \
            {input.merged_IS_spikeIn} {input.readStats_plot} {input.ISreads_plot} {input.spikeInIS_plot} \
            {input.IS_plot_per_sample} {input.IS_plot_combine_replicates} {input.IS} {input.IS_spikeIn} {input.cov} {output}
        '''

rule site__index_rmd:
    input:
        expand(SITE + '{pr}.Rmd', pr=RUNS) ## here {run} is platform--run
    params:
        index = CONFIG + 'index.Rmd',
        datasets = CONFIG + 'AAV_integration_data_analyses_overview.xlsx'
    output:
        SITE + 'index.Rmd'
    shell:
        'python bin/mk_index.py {params.index} {params.datasets} {output} {input}'


rule site__deploy_site_index:
    input:
        site = SITE + 'index.Rmd',
    log:
        'logs/deploy_site_index.log'
    output:
        touch(INTs + 'flags/index_up.' + SITE_NAME)
    container:
        'docker://samesense/tidyverse-r:5'
    shell:
        'Rscript bin/deploy_index.R {SITE} {SITE_NAME} > {log}'

rule site__target:
    input:
        INTs + 'flags/rsconnect',
        INTs + 'flags/index_up.' + SITE_NAME # deploy_site_index
