rule annofinal__IS_bed:
    input:
        INTs + "IS_results/IS_events_collapse_comb_rmcollisions/{platform}/{run}__IS.csv",
    output:
        INTs + "IS_results/IS_events_collapse_comb_rmcollisions_bed/{platform}/{run}__IS.bed",
    shell:
        '''
        if echo {input} | grep -q '230615_M05240_0126_000000000-KNTM3'; then
            echo 'not S5 - S12 type of samples, skip!!!'
        else
            csvcut --maxfieldsize=9999999 -c3,4 {input} | sed 1d | awk -F"," '{{print $1"\t"$2-1"\t"$2}}' | \
            sort -k1,1 -k2,2n | uniq | \
            sed '1ichrom\tstart\tend' > {output}
        fi
        '''

rule annofinal__closest_host_TSS:
    input:
        IS_bed = INTs + "IS_results/IS_events_collapse_comb_rmcollisions_bed/{platform}/{run}__IS.bed",
    output:
        INTs + "IS_results_final/IS_events_closest_TSS/{platform}/{run}__IS.tsv",
    params:
        host_tss_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/ncbiRefSeq_TSS.bed"
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
            # the nearest TSS could be different from the "annotated" gene, even when the IS is within a GeneID
            bedtools closest -a {input.IS_bed} -b {params.host_tss_bed} -d -D ref -t first > {output}
        '''

rule annofinal__closest_ortholog_TSS:
    input:
        IS_bed = INTs + "IS_results/IS_events_collapse_comb_rmcollisions_bed/{platform}/{run}__IS.bed",
    output:
        INTs + "IS_results_final/IS_events_closest_ortholog_TSS/{platform}/{run}__IS.tsv",
    params:
        host_ref = lambda wildcards: RUN_HOST[wildcards.run],
        human_orth_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/human_xenoRefGene_TSS.bed"
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
            if [[ {params.host_ref} != "hg38" ]]; then
                # the nearest TSS could be different from the "annotated" gene, even when the IS is within a GeneID
                bedtools closest -a {input.IS_bed} -b {params.human_orth_bed} -d -D ref -t first > {output}
            else
                touch {output}
            fi
        '''

rule annofinal__host_in_out_gene:
    input:
        IS_bed = INTs + "IS_results/IS_events_collapse_comb_rmcollisions_bed/{platform}/{run}__IS.bed",
    output:
        host_gene_InOut = INTs + "IS_results_final/IS_InOut_hostGene/{platform}/{run}__IS.tsv",
    params:
        host_tr_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/ncbiRefSeq_tr.bed",
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
        intersectBed -a {input.IS_bed} -b {params.host_tr_bed} -wao | cut -f1-3,10 | sort| uniq > {output.host_gene_InOut}
        '''


rule annofinal__ortholog_in_out_gene:
    input:
        IS_bed = INTs + "IS_results/IS_events_collapse_comb_rmcollisions_bed/{platform}/{run}__IS.bed",
    output:
        human_orth_gene_InOut = INTs + "IS_results_final/IS_InOut_humanorthGene/{platform}/{run}__IS.tsv",
    params:
        host_ref = lambda wildcards: RUN_HOST[wildcards.run],
        human_orth_tr_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/human_xenoRefGene_tr.bed"
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
        if [[ {params.host_ref} != "hg38" ]]; then
            intersectBed -a {input.IS_bed} -b {params.human_orth_tr_bed} -wao | cut -f1-3,10 | sort| uniq > {output.human_orth_gene_InOut}
        else
            touch {output}
        fi
        '''

rule annofinal__in_out_features:
    input:
        IS_bed = INTs + "IS_results/IS_events_collapse_comb_rmcollisions_bed/{platform}/{run}__IS.bed",

    output:
        liver_enhancer_InOut = INTs + "IS_results_final/IS_liver_enhancer_InOut/{platform}/{run}__IS.tsv",
        CpG_InOut = INTs + "IS_results_final/IS_CpG_InOut/{platform}/{run}__IS.tsv",
        rmsk_anno = INTs + "IS_results_final/IS_rmsk_anno/{platform}/{run}__IS.tsv"
    params:
        liver_enhancer_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/liver_enhancer.bed",
        CpG_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/CpG.bed",
        rmsk_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/rmsk.bed",
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
        intersectBed -a {input.IS_bed} -b {params.liver_enhancer_bed} -wao | cut -f1-3,7 | sort| uniq > {output.liver_enhancer_InOut}
        intersectBed -a {input.IS_bed} -b {params.CpG_bed} -wao | cut -f1-3,7 | sort| uniq > {output.CpG_InOut}
        intersectBed -a {input.IS_bed} -b {params.rmsk_bed} -wao | cut -f1-3,7 | \
            awk '{{ key = $1"_"$2"_"$3; if(key in combined) {{combined[key] = combined[key] "/" $4}} else {{combined[key] = $0}} }} \
                END {{ for(entry in combined) {{print combined[entry]}} }}' > {output.rmsk_anno}
        '''

rule annofinal__add_annotations:
    input:
        IS_events = INTs + "IS_results/IS_events_collapse_comb_rmcollisions/{platform}/{run}__IS.csv",
        host_tss_bed = INTs + "IS_results_final/IS_events_closest_TSS/{platform}/{run}__IS.tsv",
        human_orth_bed = INTs + "IS_results_final/IS_events_closest_ortholog_TSS/{platform}/{run}__IS.tsv",
        host_gene_InOut = INTs + "IS_results_final/IS_InOut_hostGene/{platform}/{run}__IS.tsv",
        human_orth_gene_InOut = INTs + "IS_results_final/IS_InOut_humanorthGene/{platform}/{run}__IS.tsv",
        liver_enhancer_InOut = INTs + "IS_results_final/IS_liver_enhancer_InOut/{platform}/{run}__IS.tsv",
        CpG_InOut = INTs + "IS_results_final/IS_CpG_InOut/{platform}/{run}__IS.tsv",
        rmsk_anno = INTs + "IS_results_final/IS_rmsk_anno/{platform}/{run}__IS.tsv"
    output:
        INTs + "IS_results_final/IS_events_anno/{platform}/{run}__IS_final.anno.csv"
    params:
        host_ref = lambda wildcards: RUN_HOST[wildcards.run]
    shell:
        '''
        if [[ {params.host_ref} != "hg38" ]]; then
            python bin/IS_anno_final.py {input.IS_events} {input.host_tss_bed} {input.human_orth_bed} {input.host_gene_InOut} {input.human_orth_gene_InOut} \
            {input.liver_enhancer_InOut} {input.CpG_InOut} {input.rmsk_anno} {output}
        else
            python bin/IS_anno_final_human.py {input.IS_events} {input.host_tss_bed} {input.host_gene_InOut} \
            {input.liver_enhancer_InOut} {input.CpG_InOut} {input.rmsk_anno} {output}
        fi
        '''


# for middle step test

def annofinal_collect_outputs_byRun(wc):
    ls = []
    for each in RUNS:
        platform, run = each.split("--")
        if run != '230615_M05240_0126_000000000-KNTM3':
            ls.append(INTs + f"IS_results_final/IS_events_anno/{platform}/{run}__IS_final.anno.csv"),
    assert len(ls) >= 1
    return ls

rule annofinal__target:
    input:
        annofinal_collect_outputs_byRun