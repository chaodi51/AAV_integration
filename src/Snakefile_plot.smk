## make plots
rule plot__read_stats:
    input:
        INTs + "{run}__mapped_table.tsv",
    output:
        INTs + "plots/{run}__readStats.html",
    shell:
        '''
        mkdir -p {INTs}plots/plot__read_stats/{wildcards.run}/
        cp bin/plot_readstats.rmd {INTs}plots/plot__read_stats/{wildcards.run}/plot_readstats.rmd

        Rscript -e '
            setwd("{INTs}plots/plot__read_stats/{wildcards.run}/");
            rmarkdown::render("{INTs}plots/plot__read_stats/{wildcards.run}/plot_readstats.rmd", \
            output_format = "html_document", params = list(map_tab = "{input}"), \
            output_file = "{output}")'
            
        rm -r {INTs}plots/plot__read_stats/{wildcards.run}/
        '''


rule plot__spikeInIS:
    input:
        map_tab = INTs + "{run}__mapped_table.tsv",
        IS_tab = INTs + "spikeIn_IS_results/spikeIn_IS_events_comb/{run}__spikeIn_IS.csv",
    output:
        INTs + "plots/{run}__spikeInIS.html",
    shell:
        '''
        mkdir -p {INTs}plots/plot__spikeInIS/{wildcards.run}/
        cp bin/plot_spikeInIS.rmd {INTs}plots/plot__spikeInIS/{wildcards.run}/plot_spikeInIS.rmd

        Rscript -e '
            setwd("{INTs}plots/plot__spikeInIS/{wildcards.run}/");
            rmarkdown::render("{INTs}plots/plot__spikeInIS/{wildcards.run}/plot_spikeInIS.rmd", \
            output_format = "html_document", params = list(map_tab = "{input.map_tab}", IS_tab = "{input.IS_tab}" ), \
            output_file = "{output}")'
            
        rm -r {INTs}plots/plot__spikeInIS/{wildcards.run}/
        '''

rule plot__IS:
    input:
        map_tab = INTs + "{run}__mapped_table.tsv",
        IS_tab = INTs + "IS_results_final/IS_events_anno/{run}__IS_final.anno.csv",
        vector_anno = DATA + 'refs/vector/anno/CAG_EGFP_ITR_corrected.bed',
        master_mapping = MAP_FILE
    output:
        per_sample = INTs + "plots/{run}__IS_per_sample.html",
        combine_replicates = INTs + "plots/{run}__IS_combine_replicates.html",
    params:
        host_ref = lambda wildcards: RUN_HOST[wildcards.run.lstrip('shortRead/')],
        n_rep = lambda wildcards: len(RUN_REP[wildcards.run.lstrip('shortRead/')])
    shell:
        '''
            ## for data does not have replicates (like vector genome run), make fake combine reps report html to run through the pipeline 
            mkdir -p {INTs}plots/plot__IS/{wildcards.run}/
            mkdir -p {INTs}plots/plot__IS_comb/{wildcards.run}/
            if [[ {params.host_ref} != "hg38" ]]; then
                cp bin/plotIS.rmd {INTs}plots/plot__IS/{wildcards.run}/plotIS.rmd
                Rscript -e '
                    setwd("{INTs}plots/plot__IS/{wildcards.run}/");
                    rmarkdown::render("{INTs}plots/plot__IS/{wildcards.run}/plotIS.rmd", \
                    output_format = "html_document", params = list(map_tab = "{input.map_tab}", IS_tab = "{input.IS_tab}", anno = "{input.vector_anno}" ), \
                    output_file = "{output.per_sample}")'
                
                if [[ {params.n_rep} -gt 1 ]]; then
                    cp bin/plotIS_combine_reps.rmd {INTs}plots/plot__IS_comb/{wildcards.run}/plotIS_combine_reps.rmd
                    Rscript -e '
                        setwd("{INTs}plots/plot__IS_comb/{wildcards.run}/");
                        rmarkdown::render("{INTs}plots/plot__IS_comb/{wildcards.run}/plotIS_combine_reps.rmd", \
                        output_format = "html_document", params = list(map_tab = "{input.map_tab}", IS_tab = "{input.IS_tab}", anno = "{input.vector_anno}", master_mapping = "{input.master_mapping}" ), \
                        output_file = "{output.combine_replicates}")'
                else
                    cp {output.per_sample}  {output.combine_replicates}
                fi

            else
                cp bin/plotIS_human.rmd {INTs}plots/plot__IS/{wildcards.run}/plotIS_human.rmd
                Rscript -e '
                    setwd("{INTs}plots/plot__IS/{wildcards.run}/");
                    rmarkdown::render("{INTs}plots/plot__IS/{wildcards.run}/plotIS_human.rmd", \
                    output_format = "html_document", params = list(map_tab = "{input.map_tab}", IS_tab = "{input.IS_tab}", anno = "{input.vector_anno}" ), \
                    output_file = "{output.per_sample}")'       
                cp {output.per_sample}  {output.combine_replicates}           
            fi
            rm -r {INTs}plots/plot__IS/{wildcards.run}/
            rm -r {INTs}plots/plot__IS_comb/{wildcards.run}/
        '''


def map_plots_output(wc):
    ls = []
    for each in RUNS:
        platform, runid = each.split("--")
        run = platform + '/' + runid
        ls.append(INTs + f"plots/{run}__readStats.html")
        ls.append(INTs + f"plots/{run}__IS_per_sample.html")
        ls.append(INTs + f"plots/{run}__IS_combine_replicates.html")
        ls.append(INTs + f"plots/{run}__spikeInIS.html")
    assert len(ls) >= 1
    return ls

rule plot__target:
    input:
        map_plots_output
