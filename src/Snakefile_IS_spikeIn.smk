rule IS_spikeIn__split_reads:
    input:
        DATA + "processed/spikeIn_plusITR_bams_rmdup/{platform}/{run}/{sample}__spikeIn_plusITR__rmdup.sorted.bam"
    output:
        sam = DATA + "processed/spikeIn_plusITR_bams_split_reads/{platform}/{run}/{sample}__spikeIn_plusITR.splitreads.sam",
    params:
        tmp = DATA + "processed/tmp/{sample}__spikeIn_plusITR.sam",
    container:
        "docker://mgibio/alignment_helper-cwl:2.2.1"
    threads: 12
    shell:
        '''
            mkdir -p "../data/processed/tmp/"
            touch {params.tmp}        
            samtools sort -n -@ {threads} {input} | samtools view -h | samblaster -e --ignoreUnmated -s {output.sam} \
            --maxSplitCount 2 --minIndelSize 30 --maxUnmappedBases 30 --minNonOverlap 20 -o {params.tmp} 
            rm -f {params.tmp}
        '''

rule IS_spikeIn__split_reads_parts:
    input:
        sam = DATA + "processed/spikeIn_plusITR_bams_split_reads/{platform}/{run}/{sample}__spikeIn_plusITR.splitreads.sam",
    output:
        pri_sam = DATA + "processed/spikeIn_plusITR_bams_split_reads_pri/{platform}/{run}/{sample}__spikeIn_plusITR.splitreads_pri.sam",
        # suppl_sam = DATA + "processed/spikeIn_plusITR_bams_split_reads_suppl/{platform}/{run}/{sample}__spikeIn_plusITR.splitreads_suppl.sam",
    container:
        "docker://samesense/bwa-samtools:v1"
    threads: 12
    shell:
        '''
            samtools view -h -F 2048 {input} > {output.pri_sam}
        '''

rule IS_spikeIn__cigar_parse:
	input:
		pri_sam = DATA + "processed/spikeIn_plusITR_bams_split_reads_pri/{platform}/{run}/{sample}__spikeIn_plusITR.splitreads_pri.sam",
	output:
		IS_tab = INTs + "spikeIn_IS_results/spikeIn_IS_reads/{platform}/{run}/{sample}__IS.csv",
		chimera_tab = INTs + "spikeIn_chimera_results/spikeIn_chimera_reads/{platform}/{run}/{sample}__chimera.csv",
	shell:
		'''
			python bin/parse_cigar.py -i {input} -o {output.IS_tab} -c {output.chimera_tab}
		'''

# format to IS events with sonication abundance
rule IS_spikeIn__format_IS:
    input:
        IS = INTs + "spikeIn_IS_results/spikeIn_IS_reads/{platform}/{run}/{sample}__IS.csv",
        bam = DATA + "processed/spikeIn_plusITR_bams_rmdup/{platform}/{run}/{sample}__spikeIn_plusITR__rmdup.sorted.bam"
    output:
        INTs + "spikeIn_IS_results/spikeIn_IS_events/{platform}/{run}/{sample}__IS.csv",
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
            python bin/sonication_sites.py -i {input.IS} -b {input.bam} -o {output} -m TES

        else
            python bin/sonication_sites.py -i {input.IS} -b {input.bam} -o {output} -m slimPCR
        fi
        '''
        
# for middle step test
def collect_outputs(wc):
	ls = []
	for each in SAMPLE_REFS:
		platform, run, sample = each.split("/")
		for ref in SAMPLE_REFS[each]:
			ls.append(INTs + f"spikeIn_IS_results/spikeIn_IS_events/{platform}/{run}/{sample}__IS.csv")
			ls.append(INTs + f"spikeIn_chimera_results/spikeIn_chimera_reads/{platform}/{run}/{sample}__chimera.csv")
	assert len(ls) >= 1
	return ls

rule IS_spikeIn__target:
    input:
        collect_outputs