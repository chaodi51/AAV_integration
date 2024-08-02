
rule ref__ref_idx:
    input:
        DATA + "refs/vector/sequence/{ref}.fa"
    output:
        DATA + "refs/vector/sequence/{ref}.fa.ann" # any of the 5 index files
    container:
        "docker://samesense/bwa-samtools:v1"    
    shell:
        'bwa index {input}'


rule ref__host_idx:
    input:
        DATA + "refs/host/{host_ref}/sequence/genome.fa"
    output:
        DATA + "refs/host/{host_ref}/sequence/genome.fa.ann"
    container:
        "docker://samesense/bwa-samtools:v1"
    shell:
        'bwa index {input}'

rule ref__combine_genome:
    input:
        vector = DATA + "refs/vector/sequence/{ref}.fa",
        host = DATA + "refs/host/{host_ref}/sequence/genome.fa"
    output:
        mixed = DATA + "refs/mixed/{host_ref}/sequence/host__{ref}.fa"
    shell:
        'cat {input.host} {input.vector} > {output.mixed}'

rule ref__mixed_idx:
    input:
        DATA + "refs/mixed/{host_ref}/sequence/host__{ref}.fa"
    output:
        DATA + "refs/mixed/{host_ref}/sequence/host__{ref}.fa.ann"
    container:
        "docker://samesense/bwa-samtools:v1"    
    shell:
        'bwa index {input}'

rule ref__spikeIn_ref:
    input:
        DATA + "refs/spikeIn/sequence/aav-is-standards_woITR.fa"
    output:
        DATA + "refs/spikeIn/sequence/aav-is-standards_woITR.fa.ann"
    container:
        "docker://samesense/bwa-samtools:v1"
    shell:
        'bwa index {input}'

rule ref__spikeIn_plusITR_ref:
    input:
        DATA + "refs/spikeIn/sequence/aav-is-standards_plusITR.fa"
    output:
        DATA + "refs/spikeIn/sequence/aav-is-standards_plusITR.fa.ann"
    container:
        "docker://samesense/bwa-samtools:v1"    
    shell:
        'bwa index {input}'

rule ref__mk_anno:
    input:
        vector_anno = DATA + "refs/vector/anno/{ref}.gtf",
        host_anno = DATA + "refs/host/{host_ref}/anno/genes.gtf"
    output:
        DATA + 'refs/mixed/{host_ref}/anno/host__{ref}.gtf'
    shell:
        'cat {input.host_anno} {input.vector_anno} > {output}'

rule ref__TSS_bed:
    input:
        DATA + "refs/host/{host_ref}/anno/ncbiRefSeq.txt"
    output:
        DATA + "refs/host/{host_ref}/anno/ncbiRefSeq_TSS.bed"
    shell:
        '''
        awk '$4=="+" {{print $3"\t"$5"\t"$5+1"\t"$13";"$2";"$3";"$5";"$6";"$4"\t.\t"$4}} \
            $4=="-" {{print $3"\t"$6-1"\t"$6"\t"$13";"$2";"$3";"$5";"$6";"$4"\t.\t"$4}}' {input} | \
            sort -k1,1 -k2,2n > {output}
        '''
rule ref__tr_bed:
    input:
        DATA + "refs/host/{host_ref}/anno/ncbiRefSeq.txt"
    output:
        DATA + "refs/host/{host_ref}/anno/ncbiRefSeq_tr.bed"
    shell:
        '''
        awk '{{print $3"\t"$5"\t"$6"\t"$13";"$2"\t.\t"$4}}' {input} | \
            sort -k1,1 -k2,2n > {output}
        '''

rule ref__mouse_human_xenoRef:
    input:
        mouse_xenoref = DATA + "refs/host/mm39/anno/xenoRefGene.txt",
        human_ncbirefseq = DATA + "refs/host/human/anno/ncbiRefSeq.txt"
    output:
        DATA + "refs/host/mm39/anno/human_xenoRefGene.txt"
    shell:
        '''
        # extract human genes from xenoRef
        python bin/extract_ids.py {input.human_ncbirefseq} {input.mouse_xenoref} {output}
        '''

rule ref__mouse_human_xenoRef_TSS_bed:
    input:
        DATA + "refs/host/mm39/anno/human_xenoRefGene.txt"
    output:
        DATA + "refs/host/mm39/anno/human_xenoRefGene_TSS.bed"
    shell:
        '''
        awk '$4=="+" {{print $3"\t"$5"\t"$5+1"\t"$13";"$2";"$3";"$5";"$6";"$4"\t.\t"$4}} \
            $4=="-" {{print $3"\t"$6-1"\t"$6"\t"$13";"$2";"$3";"$5";"$6";"$4"\t.\t"$4}}' {input} | \
            sort -k1,1 -k2,2n > {output}
        '''

rule ref__mouse_human_xenoRef_tr_bed:
    input:
        DATA + "refs/host/mm39/anno/human_xenoRefGene.txt"
    output:
        DATA + "refs/host/mm39/anno/human_xenoRefGene_tr.bed"
    shell:
        '''
        awk '{{print $3"\t"$5"\t"$6"\t"$13";"$2"\t.\t"$4}}' {input} | \
            sort -k1,1 -k2,2n > {output}
            '''

rule ref__mk_regions:
    params:
        bed = DATA + "refs/vector/master_anno.bed", # cat all vector.bed
    output:
        regions = INTs + "regions.standardized.tsv",
    run:
        import pandas as pd
        cols = ["chrom", "st", "end", "region", "val", "strand"]
        df = pd.read_csv(params.bed, sep="\t", names=cols).rename(
            columns={"st": "st_pre"}
        )
        crit = df.apply(
            lambda row: (row["end"] - row["st_pre"]) > 50
            and len(str(row["region"])) < 50,
            axis=1,
        )
        df = df[crit]
        df.loc[:, "st"] = df["st_pre"] + 1

        refs = set(df["chrom"])
        dummy_df = []
        for ref in REFS:
            if not ref in refs:
                row = {
                    "chrom": ref,
                    "st": 1,
                    "end": 1,
                    "region": "dummy",
                    "val": 0,
                    "strand": ".",
                }
                dummy_df.append(row)

        df = pd.concat([df, pd.DataFrame(dummy_df, columns=cols)])
        cols2 = ["chrom", "st", "end", "region", "val", "strand"]
        df[cols2].to_csv(output.regions, sep="\t", index=False)


def target_input_ref(wc):
    ls = []
    ls.append(DATA + f"refs/spikeIn/sequence/aav-is-standards.fa.ann")
    ls.append(DATA + f"refs/spikeIn/sequence/aav-is-standards_plusITR.fa.ann")
    ls.append(DATA + f"refs/host/mm39/anno/human_xenoRefGene_TSS.bed")
    ls.append(DATA + f"refs/host/mm39/anno/human_xenoRefGene_tr.bed")
    ls.append(INTs + f"regions.standardized.tsv")

    for host_ref in HOST_REF:
        ls.append(DATA + f"refs/host/{host_ref}/anno/ncbiRefSeq_TSS.bed")
        ls.append(DATA + f"refs/host/{host_ref}/anno/ncbiRefSeq_tr.bed")
        ls.append(DATA + f"refs/host/{host_ref}/sequence/genome.fa.ann")
        for ref in REFS:
            ls.append(DATA + f"refs/mixed/{host_ref}/sequence/host__{ref}.fa")
            ls.append(DATA + f"refs/mixed/{host_ref}/sequence/host__{ref}.fa.ann")
            ls.append(DATA + f'refs/mixed/{host_ref}/anno/host__{ref}.gtf')
            ls.append(DATA + f"refs/vector/sequence/{ref}.fa.ann")

    assert len(ls) >= 1
    return ls

# for middle processing test
rule ref__target:
    input:
        target_input_ref
