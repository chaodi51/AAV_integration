## count total reads
rule count__tot_init_reads:
    input:
        nolinker_r1 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R1.fastq.gz",
    output:
        INTs + "tot_reads/{platform}/{run}/{sample}.tot.count",
    shell:
        "gunzip -c {input} | grep '^@' | awk '(NF>1){{print}}' - "
        "| cut -f 1 -d ' ' | cut -f 2 -d '@' | wc -l > {output}"

## count reads mapped to spike-in 
rule count__tot_reads_spikeIn:
    input:
        bam = DATA + "processed/spikeIn_plusITR_bams_rmdup/{platform}/{run}/{sample}__spikeIn_plusITR__rmdup.sorted.bam"
    output:
        INTs + "bamwork/mapped_count_spikeIn/{platform}/{run}/{sample}__{ref}.tot.count",
    shell:
        '''
        samtools flagstats {input.bam} |grep "primary mapped" | \
        cut -d " " -f1 > {output}
        '''

rule count__each_reads_spikeIn:
    input:
        bam = DATA + "processed/spikeIn_plusITR_bams_rmdup/{platform}/{run}/{sample}__spikeIn_plusITR__rmdup.sorted.bam"
    output:
        INTs + "bamwork/mapped_count_spikeIn/{platform}/{run}/{sample}__{ref}.each.count",
    shell:
        '''
        samtools view -b -F 2048 {input.bam} | samtools idxstats - | head -n -1 | cut -f3 | paste -s -d',' - > {output}
        '''

## count reads mapped to vector 
rule count__reads_vector:
    input:
        bam = DATA + "processed/vector_bams/{platform}/{run}/{sample}__{ref}.sorted.bam"
    output:
        INTs + "bamwork/mapped_count_vector/{platform}/{run}/{sample}__{ref}.tot.count",
    shell:
        '''
        samtools flagstats {input.bam} |grep "primary mapped" | \
        cut -d " " -f1 > {output}
        '''

## count reads mapped to mixed genome
rule count__reads_mixed:
    input:
        bam = DATA + "processed/mixed_bams/{platform}/{run}/{sample}__{ref}__mixed.sorted.bam"
    output:
        INTs + "bamwork/mapped_count_mixed/{platform}/{run}/{sample}__{ref}.tot.count",
    shell:
        '''
        samtools flagstats {input.bam} |grep "primary mapped" | \
        cut -d " " -f1 > {output}
        '''

## count reads mapped to mixed genome, after remove duplicates
rule count__reads_rmdup:
    input:
        bam = DATA + "processed/mixed_bams_rmdup/{platform}/{run}/{sample}__{ref}__mixed__rmdup.sorted.bam"
    output:
        INTs + "bamwork/mapped_count_rmdup/{platform}/{run}/{sample}__{ref}.tot.count",
    shell:
        '''
        samtools flagstats {input.bam} |grep "primary mapped" | \
        cut -d " " -f1 > {output}
        '''

## count reads mapped to mixed genome, after remove duplicates and keep properly paired reads
rule count__reads_mixed_filtered:
    input:
        bam = DATA + "processed/mixed_bams_paired/{platform}/{run}/{sample}__{ref}__mixed_paired.sorted.bam"
    output:
        category_counts = INTs + "bamwork/mapped_count_rmdup/{platform}/{run}/{sample}__{ref}__mixed.category.count",
    shell:
        '''
        python bin/count_reads.py {input.bam} {output.category_counts}
        '''

def get_read_int(afile):
    with open(afile) as f:
        base = f.readline().strip()
    return base

rule count__read_df:
    input:
        tot_counts = INTs + "tot_reads/{platform}/{run}/{sample}.tot.count",
        # mapped_spikeIn_tot = INTs + "bamwork/mapped_count_spikeIn/{platform}/{run}/{sample}__{ref}.tot.count",
        mapped_spikeIn_each = INTs + "bamwork/mapped_count_spikeIn/{platform}/{run}/{sample}__{ref}.each.count",
        mapped_vector_tot = INTs + "bamwork/mapped_count_vector/{platform}/{run}/{sample}__{ref}.tot.count",
        mapped_mixed_tot = INTs + "bamwork/mapped_count_mixed/{platform}/{run}/{sample}__{ref}.tot.count",
        mapped_mixed_rmdup = INTs + "bamwork/mapped_count_rmdup/{platform}/{run}/{sample}__{ref}.tot.count",
        category_counts = INTs + "bamwork/mapped_count_rmdup/{platform}/{run}/{sample}__{ref}__mixed.category.count",
    output:
        o = INTs + "mapped_df/{platform}/{run}/{sample}__{ref}__mixed.bwa",
    run:
        import pandas as pd
        (standard1, standard2, standard3, standard4, standard5, standard6, chrV) = get_read_int(input.mapped_spikeIn_each).split(',')
        (mapped_mixed_paired, vector_reads, host_reads, chimeric_reads, vector_vector_reads, 
            host_host_reads, vector_host_reads) = get_read_int(input.category_counts).split(',')
        data = {
            "run": [wildcards.run],
            "Sample": [wildcards.sample],
            "ref": [wildcards.ref],
            "tot_input": [int(get_read_int(input.tot_counts)) * 2],
            # "tot_spikeIn": [int(get_read_int(input.mapped_spikeIn_tot))],
            "standard1": [int(standard1)],
            "standard2": [int(standard2)],
            "standard3": [int(standard3)],
            "standard4": [int(standard4)],
            "standard5": [int(standard5)],
            "standard6": [int(standard6)],
            "vectorGenome": [int(get_read_int(input.mapped_vector_tot))],
            "mapped_frac": [int(get_read_int(input.mapped_vector_tot)) / (int(get_read_int(input.tot_counts)) * 2)],
            "mixedGenome": [int(get_read_int(input.mapped_mixed_tot))],
            "mixedGenome_rmdup": [int(get_read_int(input.mapped_mixed_rmdup))],
            "mixedGenome_rmdup_paired": [int(mapped_mixed_paired)],
            "vector_containing": [int(vector_reads)],
            "host_containing": [int(host_reads)],
            "chimeric": [int(chimeric_reads)],
            "vector_vector": [int(vector_vector_reads)],
            "host_host": [int(host_host_reads)],
            "vector_host": [int(vector_host_reads)]
        }

        # Create a DataFrame
        df = pd.DataFrame(data)
        # pd.set_option('display.max_colwidth', 30)
        # Save the DataFrame to a CSV file
        df.to_csv(output.o, sep="\t", index=False)


def mk_cat_map_table_input(wc):
    ls = []
    for sample in SAMPLES:
        platform, run, name = sample.split("/")
        key = f"{platform}/{run}"
        if key == wc.run:
            for ref in SAMPLE_REFS[sample]:
                ls.append(INTs + f"mapped_df/{sample}__{ref}__mixed.bwa")
    assert len(ls) >= 1
    return ls

rule count__map_table:
    input:
        mk_cat_map_table_input,
    output:
        INTs + "{run}__mapped_table.tsv",
    conda:
        REQS + "cat.conda.env.yaml"
    shell:
        "python-cath {input} {output}"



# for middle step test
def mk_mapped_table(wc):
    ls = []
    for each in RUNS:
        platform, run = each.split("--")
        ls.append(INTs + f"{platform}/{run}__mapped_table.tsv")

    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        for ref in SAMPLE_REFS[each]:
            ls.append(INTs + f"bamwork/mapped_count_spikeIn/{platform}/{run}/{sample}__{ref}.tot.count")
            ls.append(INTs + f"bamwork/mapped_count_spikeIn/{platform}/{run}/{sample}__{ref}.each.count")

    assert len(ls) >= 1
    return ls

rule count__target:
    input:
        mk_mapped_table,