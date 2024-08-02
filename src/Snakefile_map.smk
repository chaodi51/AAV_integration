# filter R1 contain linker (ACCTAACTGCTGTGCCACTG) for PCR-based methods
# Azemta linker (AGTCCCTTAAGCGGAG), their PCR direction is opposite to us, already switched R1 and R2 beforehand.
# many linkers attached to ITR sequences (due to linear AAV genome), failed to amplify host genome
rule map__filter_R1_linker:
    input:
        r1 = TRIM_FQ_DIR + "{platform}/{run}/{sample}_R1.fastq.gz",
        r2 = TRIM_FQ_DIR + "{platform}/{run}/{sample}_R2.fastq.gz",
    output:
        nolinker_r1 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R1.fastq.gz",
        nolinker_r2 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R2.fastq.gz",
    params:
        all_read = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_all.readid",
        linker_read = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_linker.readid",
        nolinker_read = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_nolinker.readid",
        lib_type = lambda wildcards: RUN_LIB_TYPE[wildcards.run]
    shell:
        '''
#         if echo {input} | grep -q '230818_VH00163_52_AACWHTHM5\|230925_M07542_0048_000000000-L6GCD\|231114_VH00163_66_AAF55KMM5\|\
# 231215_VH00163_75_AAF2M5TM5\|240112_VH00163_77_AAF5V37M5\|230816_M71053_29_000000000-LBGW4\|\
# 240209_VH00163_85_AAF5WC2M5\|240214_VH00163_87_AAF5W7NM5\|230817_LH00401_55_223JYKLT4_TES\|240329_VH00163_95_AAFCNC7M5\|\
# 240405_M05240_0149_000000000-L9FLR\|240328_VH00163_94_AAFJHJGM5';

        if [[ {params.lib_type} == "TES" ]];
        then
            echo -e "No R2 specific linker, might be TES method!\n"
            cp {input.r1} {output.nolinker_r1}
            cp {input.r2} {output.nolinker_r2} 
        else
            ## TODO: add info to the master_mapping.tsv file and use as parameters here: machine id @XXXX, linker sequence etc.
            # LC_ALL=C # sets the character encoding to the ASCII character set
            zcat {input.r1} | grep "^@M05240\|^@M07542\|^@VH00163\|^@A01959\|^@LH00401" | cut -d " " -f1 > {params.all_read}
            zcat {input.r1} | grep -B 1 "ACCTAACTGCTGTGCCACTG\|AGTCCCTTAAGCGGAG" | grep "^@M05240\|^@M07542\|^@VH00163\|^@A01959\|^@LH00401" | cut -d " " -f1 > {params.linker_read}
            python bin/extract_readid.py {params.all_read} {params.linker_read} {params.nolinker_read}
            zcat {input.r1} | grep -A 3 -Ff {params.nolinker_read} | grep -v '^--$' | pigz -f -c > {output.nolinker_r1}
            zcat {input.r2} | grep -A 3 -Ff {params.nolinker_read} | grep -v '^--$' | pigz -f -c > {output.nolinker_r2}
        fi
        '''

# map reads to spike-in
rule map__bwa_spikeIn:
    input:
        nolinker_r1 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R1.fastq.gz",
        nolinker_r2 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R2.fastq.gz",
        ref = DATA + "refs/spikeIn/sequence/aav-is-standards_woITR.fa"
    output:
        DATA + "processed/spikeIn_bams/{platform}/{run}/{sample}__spikeIn.sorted.bam"
    container:
        "docker://samesense/bwa-samtools:v1"
    threads: 8
    shell:
        '''
        bwa mem -t {threads} {input.ref} {input.nolinker_r1} {input.nolinker_r2} | \
        samtools view -b -F 772 | \
        samtools sort - -@ {threads} -o {output} && samtools index {output}
        '''

rule map__bwa_spikeIn_plusITR:
    input:
        nolinker_r1 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R1.fastq.gz",
        nolinker_r2 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R2.fastq.gz",
        ref = DATA + "refs/spikeIn/sequence/aav-is-standards_plusITR.fa"
    output:
        DATA + "processed/spikeIn_plusITR_bams/{platform}/{run}/{sample}__spikeIn_plusITR.sorted.bam"
    container:
        "docker://samesense/bwa-samtools:v1"
    threads: 8
    shell:
        '''
        bwa mem -t {threads} {input.ref} {input.nolinker_r1} {input.nolinker_r2} | \
        samtools view -b -F 772 | \
        samtools sort - -@ {threads} -o {output} && samtools index {output}
        '''
        
rule map__spikeIn_remove_dup:
    input:
        DATA + "processed/spikeIn_plusITR_bams/{platform}/{run}/{sample}__spikeIn_plusITR.sorted.bam"
    output:
        DATA + "processed/spikeIn_plusITR_bams_rmdup/{platform}/{run}/{sample}__spikeIn_plusITR__rmdup.sorted.bam"
    container:
        "docker://lethalfang/sambamba:1.0.0"
    threads: 8
    params:
        tmpdir = DATA + "processed/rmdup_tmp/{platform}/{run}/{sample}"
    shell:
        '''
        mkdir -p "../data/processed/rmdup_tmp/{wildcards.platform}/{wildcards.run}/{wildcards.sample}"
        sambamba markdup -r -t {threads} --tmpdir {params.tmpdir} {input} {output}
        ''' 

# map filtered fastq files to vector
rule map__bwa_vector:
    input:
        nolinker_r1 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R1.fastq.gz",
        nolinker_r2 = DATA + "processed/no_linkerR1_fastq/{platform}/{run}/{sample}_R2.fastq.gz",
        ref = DATA + "refs/vector/sequence/{ref}.fa"
    output:
        DATA + "processed/vector_bams/{platform}/{run}/{sample}__{ref}.sorted.bam"
    container:
        "docker://samesense/bwa-samtools:v1"
    threads: 8
    shell:
        '''
        # Montini lab used -k 16 -r 1 -M -v 1 -T 15, it maybe too permissive
        bwa mem -t {threads} {input.ref} {input.nolinker_r1} {input.nolinker_r2} | \
        samtools view -b -F 772 | \
        samtools sort - -@ {threads} -o {output} && samtools index {output}
        '''

rule map__bam_to_bed:
    input:
        DATA + "processed/vector_bams/{platform}/{run}/{sample}__{ref}.sorted.bam"
    output:
        bed = DATA + "processed/vector_beds/{platform}/{run}/{sample}__{ref}.sorted.bed",
        rid = DATA + "processed/vector_reads/{platform}/{run}/{sample}__{ref}.sorted.readid"
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        '''
        bedtools bamtobed -cigar -i {input} > {output.bed}
        cut -f4 {output.bed} | sed 's/\/.$//g' | sort |uniq > {output.rid}
        '''

rule map__extract_fq:
    input:
        rid = DATA + "processed/vector_reads/{platform}/{run}/{sample}__{ref}.sorted.readid",
        r1 = TRIM_FQ_DIR + "{platform}/{run}/{sample}_R1.fastq.gz",
        r2 = TRIM_FQ_DIR + "{platform}/{run}/{sample}_R2.fastq.gz"
    output:
        vector_r1 = DATA + "processed/vector_fqs/{platform}/{run}/{sample}__{ref}_R1.fastq.gz",
        vector_r2 = DATA + "processed/vector_fqs/{platform}/{run}/{sample}__{ref}_R2.fastq.gz"
    shell:
        '''
        zcat {input.r1} | grep -A 3 -Ff {input.rid} | grep -v '^--$' | pigz -f -c > {output.vector_r1}
        zcat {input.r2} | grep -A 3 -Ff {input.rid} | grep -v '^--$' | pigz -f -c > {output.vector_r2}
        '''

rule map__bwa_mixed:
    input:
        vector_r1 = DATA + "processed/vector_fqs/{platform}/{run}/{sample}__{ref}_R1.fastq.gz",
        vector_r2 = DATA + "processed/vector_fqs/{platform}/{run}/{sample}__{ref}_R2.fastq.gz",
    output:
        DATA + "processed/mixed_bams/{platform}/{run}/{sample}__{ref}__mixed.sorted.bam"
    params:
        ref = lambda wildcards: DATA + f"refs/mixed/{RUN_HOST[wildcards.run]}/sequence/host__{wildcards.ref}.fa"
    container:
        "docker://samesense/bwa-samtools:v1"
        # "docker://lethalfang/sambamba:1.0.0"
    threads: 8
    shell:
        '''
        # only filter non-primary mapping by -F 772, these are tagged as XA:, so no secondary alignment in the output for bwa
        bwa mem -t {threads} {params.ref} {input.vector_r1} {input.vector_r2} | \
        samtools view -b | \
        samtools sort - -@ {threads} -o {output} && samtools index {output}
        '''

rule map__remove_dup:
    input:
        DATA + "processed/mixed_bams/{platform}/{run}/{sample}__{ref}__mixed.sorted.bam"
    output:
        DATA + "processed/mixed_bams_rmdup/{platform}/{run}/{sample}__{ref}__mixed__rmdup.sorted.bam"
    container:
        "docker://lethalfang/sambamba:1.0.0"
    threads: 8
    params:
        tmpdir = DATA + "processed/rmdup_tmp/{platform}/{run}/{sample}__{ref}"
    shell:
        '''
        mkdir -p "../data/processed/rmdup_tmp/{wildcards.platform}/{wildcards.run}/{wildcards.sample}__{wildcards.ref}"
        sambamba markdup -r -t {threads} --tmpdir {params.tmpdir} {input} {output}
        '''   


# keep properly paired reads if both R1 and R2 map to host
rule map__filter_proper_paired:
    input:
        DATA + "processed/mixed_bams_rmdup/{platform}/{run}/{sample}__{ref}__mixed__rmdup.sorted.bam"
    output:
        DATA + "processed/mixed_bams_paired/{platform}/{run}/{sample}__{ref}__mixed_paired.sorted.bam"
    # container:
    #     "docker://samesense/bwa-samtools:v1"
    threads: 6
    shell:
        '''
        python bin/filter_host_proper_paired_reads.py -i {input} -o {output}
        '''

## filter unique reads ONLY for host, not used in IS detection, just for test !!!!
rule map__uniqe_host_align:
    input:
        DATA + "processed/mixed_bams_rmdup/{platform}/{run}/{sample}__{ref}__mixed__rmdup.sorted.bam"
    output:
        uniq_bam = DATA + "processed/mixed_bams_unique/{platform}/{run}/{sample}__{ref}__mixed__unique.sorted.bam",
        uniq_label = DATA + "processed/mixed_bams_host_unique_label/{platform}/{run}/{sample}__{ref}__mixed__host_unique_label.csv"
    shell:
        '''
        python bin/filter_and_label_unique.py -i {input} -o {output.uniq_bam} -u {output.uniq_label}
        '''    
        
# for middle step test
def collect_outputs(wc):
    ls = []
    for full_sample, ref_dict in SAMPLE_REFS.items():
        platform, run, sample = full_sample.split("/")
        # ls.append(DATA + f"processed/spikeIn_bams/{platform}/{run}/{sample}__spikeIn.sorted.bam")
        ls.append(DATA + f"processed/spikeIn_plusITR_bams/{platform}/{run}/{sample}__spikeIn_plusITR.sorted.bam")
        ls.append(DATA + f"processed/spikeIn_plusITR_bams_rmdup/{platform}/{run}/{sample}__spikeIn_plusITR__rmdup.sorted.bam")
        for ref in ref_dict:
            ls.append(DATA + f"processed/mixed_bams_host_unique_label/{platform}/{run}/{sample}__{ref}__mixed__host_unique_label.csv")
            ls.append(DATA + f"processed/mixed_bams/{platform}/{run}/{sample}__{ref}__mixed.sorted.bam")
            ls.append(DATA + f"processed/mixed_bams_paired/{platform}/{run}/{sample}__{ref}__mixed_paired.sorted.bam")

    assert len(ls) >= 1
    return ls

rule map__target:
    input:
        collect_outputs
