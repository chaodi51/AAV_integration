## merge IS within 10bp and remove collisions between un-related samples within a run

rule cleanIS__cluster:
    input:
        INTs + "IS_results/IS_events_anno/{platform}/{run}/{sample}__{ref}__IS.csv"
    output:
        INTs + "IS_results/IS_events_collapse/{platform}/{run}/{sample}__{ref}__IS.csv"
    shell:
        '''
        python bin/collapseIS.py {input} {output}
        '''

def list_IScluster_tables(wc):
    ls = []
    for sample in SAMPLES:
        platform, run, name = sample.split("/")
        key = f"{platform}/{run}"
        if key == wc.run:
            for ref in SAMPLE_REFS[sample]:
                ls.append(INTs + f"IS_results/IS_events_collapse/{sample}__{ref}__IS.csv")
    assert len(ls) >= 1
    return ls

rule cleanIS__merge_tab:
    input:
        list_IScluster_tables
    output:
        INTs + "IS_results/IS_events_collapse_comb/{run}__IS.csv"
    conda:
        REQS + "cat.conda.env.yaml"
    shell:
        "python-cath {input} {output}"

# for middle step test
def list_samples(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        for ref in SAMPLE_REFS[each]:
            ls.append(INTs + f"IS_results/IS_events_collapse/{platform}/{run}/{sample}__{ref}__IS.csv")
    assert len(ls) >= 1
    return ls



rule cleanIS__rm_collisions:
    input:
        INTs + "IS_results/IS_events_collapse_comb/{run}__IS.csv"
    output:
        INTs + "IS_results/IS_events_collapse_comb_rmcollisions/{run}__IS.csv"
    params:
        map_file = MAP_FILE,
    shell:
        '''
        if echo {input} | grep -q '230615_M05240_0126_000000000-KNTM3\|240328_VH00163_94_AAFJHJGM5\|240405_M05240_0149_000000000-L9FLR'; then
            echo 'not S5 - S12 type of samples, skip rm_collisions!!!'
            python bin/rm_collisions.py {wildcards.run} {params.map_file} {input} {output}
        else 
            python bin/rm_collisions.py {wildcards.run} {params.map_file} {input} {output}
        fi
        '''

def list_runs(wc):
    ls = []
    for each in RUNS:
        platform, run = each.split("--")
        ls.append(INTs + f"IS_results/IS_events_collapse_comb/{platform}/{run}__IS.csv"),
        if run != '230615_M05240_0126_000000000-KNTM3':
            ls.append(INTs + f"IS_results/IS_events_collapse_comb_rmcollisions/{platform}/{run}__IS.csv")
    assert len(ls) >= 1
    return ls
    
rule cleanIS__target:
    input:
        list_runs
