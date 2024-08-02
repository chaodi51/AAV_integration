rule anno__IS_bed:
    input:
        INTs + "IS_results/IS_events/{platform}/{run}/{sample}__{ref}__IS.csv",
    output:
        INTs + "IS_results/IS_events_bed/{platform}/{run}/{sample}__{ref}__IS.bed",
    shell:
        '''
        csvcut --maxfieldsize=9999999 -c3,4 {input} | sed 1d | awk -F"," '{{print $1"\t"$2-1"\t"$2}}' | \
        sort -k1,1 -k2,2n | \
        sed '1ichrom\tstart\tend' > {output}
        '''

rule anno__closest_host_TSS:
    input:
        IS_bed = INTs + "IS_results/IS_events_bed/{platform}/{run}/{sample}__{ref}__IS.bed",
    output:
        INTs + "IS_results/IS_closest_TSS/{platform}/{run}/{sample}__{ref}__IS.tsv",
    params:
        host_tss_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/ncbiRefSeq_TSS.bed"
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
            # the nearest TSS could be different from the "annotated" gene, even when the IS is within a GeneID
            bedtools closest -a {input.IS_bed} -b {params.host_tss_bed} -d -D ref -t first > {output}
        '''

rule anno__closest_ortholog_TSS:
    input:
        IS_bed = INTs + "IS_results/IS_events_bed/{platform}/{run}/{sample}__{ref}__IS.bed",
    output:
        INTs + "IS_results/IS_closest_ortholog_TSS/{platform}/{run}/{sample}__{ref}__IS.tsv",
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

rule anno__host_in_out_gene:
    input:
        IS_bed = INTs + "IS_results/IS_events_bed/{platform}/{run}/{sample}__{ref}__IS.bed",
        
    output:
        host_gene_InOut = INTs + "IS_results/IS_InOut_hostGene/{platform}/{run}/{sample}__{ref}__IS.tsv",
    params:
        host_tr_bed = lambda wildcards: DATA + f"refs/host/{RUN_HOST[wildcards.run]}/anno/ncbiRefSeq_tr.bed",
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
        intersectBed -a {input.IS_bed} -b {params.host_tr_bed} -wao | cut -f1-3,10 | sort| uniq > {output.host_gene_InOut}
        '''

rule anno__ortholog_in_out_gene:
    input:
        IS_bed = INTs + "IS_results/IS_events_bed/{platform}/{run}/{sample}__{ref}__IS.bed",
        
    output:
        human_orth_gene_InOut = INTs + "IS_results/IS_InOut_humanorthGene/{platform}/{run}/{sample}__{ref}__IS.tsv",
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


rule anno__in_out_features:
    input:
        IS_bed = INTs + "IS_results/IS_events_bed/{platform}/{run}/{sample}__{ref}__IS.bed",

    output:
        liver_enhancer_InOut = INTs + "IS_results/IS_liver_enhancer_InOut/{platform}/{run}/{sample}__{ref}__IS.tsv",
        CpG_InOut = INTs + "IS_results/IS_CpG_InOut/{platform}/{run}/{sample}__{ref}__IS.tsv",
        rmsk_anno = INTs + "IS_results/IS_rmsk_anno/{platform}/{run}/{sample}__{ref}__IS.tsv"
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

rule anno__add_annotations:
    input:
        IS_events = INTs + "IS_results/IS_events/{platform}/{run}/{sample}__{ref}__IS.csv",
        host_tss_bed = INTs + "IS_results/IS_closest_TSS/{platform}/{run}/{sample}__{ref}__IS.tsv",
        human_orth_bed = INTs + "IS_results/IS_closest_ortholog_TSS/{platform}/{run}/{sample}__{ref}__IS.tsv",
        host_gene_InOut = INTs + "IS_results/IS_InOut_hostGene/{platform}/{run}/{sample}__{ref}__IS.tsv",
        human_orth_gene_InOut = INTs + "IS_results/IS_InOut_humanorthGene/{platform}/{run}/{sample}__{ref}__IS.tsv",
        liver_enhancer_InOut = INTs + "IS_results/IS_liver_enhancer_InOut/{platform}/{run}/{sample}__{ref}__IS.tsv",
        CpG_InOut = INTs + "IS_results/IS_CpG_InOut/{platform}/{run}/{sample}__{ref}__IS.tsv",
        rmsk_anno = INTs + "IS_results/IS_rmsk_anno/{platform}/{run}/{sample}__{ref}__IS.tsv"
    output:
        INTs + "IS_results/IS_events_anno/{platform}/{run}/{sample}__{ref}__IS.csv"
    params:
        host_ref = lambda wildcards: RUN_HOST[wildcards.run]
    shell:
        '''
        if [[ {params.host_ref} != "hg38" ]]; then
            python bin/IS_anno_v1.py {input.IS_events} {input.host_tss_bed} {input.human_orth_bed} {input.host_gene_InOut} {input.human_orth_gene_InOut} \
            {input.liver_enhancer_InOut} {input.CpG_InOut} {input.rmsk_anno} {output}
        else
            python bin/IS_anno_v1_human.py {input.IS_events} {input.host_tss_bed} {input.host_gene_InOut} \
            {input.liver_enhancer_InOut} {input.CpG_InOut} {input.rmsk_anno} {output}
        fi
        '''


def list_IS_tables(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        if run == wc.run:
            for ref in SAMPLE_REFS[each]:
                ls.append(INTs + f"IS_results/IS_events_anno/{platform}/{run}/{sample}__{ref}__IS.csv")
    assert len(ls) >= 1
    return ls

rule anno__merge_tab:
    input:
        list_IS_tables
    output:
        INTs + "IS_results/IS_events_comb/{platform}/{run}__IS.anno.csv"
    conda:
        REQS + "cat.conda.env.yaml"
    shell:
        "python-cath {input} {output}"

## for spikeIn IS
rule anno__add_annotations_spikeIn:
    input:
        INTs + "spikeIn_IS_results/spikeIn_IS_events/{platform}/{run}/{sample}__IS.csv",
    output:
        INTs + "spikeIn_IS_results/spikeIn_IS_events_anno/{platform}/{run}/{sample}__IS.csv"
    shell:
        '''
            header=$(head -1 {input})
            new_header=`echo "sample",$header`
            awk -v sample={wildcards.sample} 'BEGIN{{FS=OFS=","}}; NR==1{{print "'$new_header'"}} NR>1{{print sample,$0}}' {input} > {output}
        '''

def list_IS_tables_spikeIn(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        if run == wc.run:
            ls.append(INTs + f"spikeIn_IS_results/spikeIn_IS_events_anno/{platform}/{run}/{sample}__IS.csv")
    assert len(ls) >= 1
    return ls

rule anno__merge_tab_spikeIn:
    input:
        list_IS_tables_spikeIn
    output:
        INTs + "spikeIn_IS_results/spikeIn_IS_events_comb/{platform}/{run}__spikeIn_IS.csv"
    conda:
        REQS + "cat.conda.env.yaml"
    shell:
        """
        python-cath {input} {output}
        """

# for middle step test
def collect_outputs_bySample(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        ls.append(INTs + f"spikeIn_IS_results/spikeIn_IS_events_anno/{platform}/{run}/{sample}__IS.csv")
        for ref in SAMPLE_REFS[each]:
            # ls.append(INTs + f"IS_results/IS_InOut_hostGene/{platform}/{run}/{sample}__{ref}__IS.tsv")
            # ls.append(INTs + f"IS_results/IS_InOut_humanorthGene/{platform}/{run}/{sample}__{ref}__IS.tsv")
            # ls.append(INTs + f"IS_results/IS_rmsk_anno/{platform}/{run}/{sample}__{ref}__IS.tsv")
            ls.append(INTs + f"IS_results/IS_events_anno/{platform}/{run}/{sample}__{ref}__IS.csv")
    assert len(ls) >= 1
    return ls

def collect_outputs_byRun(wc):
    ls = []
    for each in RUNS:
        platform, run = each.split("--")
        ls.append(INTs + f"IS_results/IS_events_comb/{platform}/{run}__IS.anno.csv"),
        ls.append(INTs + f"spikeIn_IS_results/spikeIn_IS_events_comb/{platform}/{run}__spikeIn_IS.csv")
    assert len(ls) >= 1
    return ls

rule anno__target:
    input:
        collect_outputs_bySample,
        collect_outputs_byRun