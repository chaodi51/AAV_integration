rule cov__dat_tmp:
    input:
        bam = DATA + "processed/vector_bams/{platform}/{run}/{sample}__{ref}.sorted.bam"
    output:
        INTs + "processed/cov_tmp/{platform}/{run}/{sample}__{ref}.cov",
    container:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        "bedtools genomecov -split -d -ibam {input.bam} | cut -f2- -| sed '1ipos\tdepth' > {output}"


rule cov__combo_depths_name:
    input:
        i = INTs + "processed/cov_tmp/{platform}/{run}/{sample}__{ref}.cov",
    output:
        o = INTs + "processed/cov_tmp/{platform}/{run}/{sample}__{ref}.finalCov.named.tsv"
    run:
        with open(input.i) as f, open(output.o, "w") as fout:
            print(f.readline().strip() + "\tplatform\trun\tsample_name", file=fout)
            for line in f:
                sample_name = wildcards.sample
                print(
                    line.strip()
                    + f"\t{wildcards.platform}\t{wildcards.run}\t{sample_name}",
                    file=fout,
                )

def mk_cov_input(wc):
    ls = []
    for each in SAMPLE_REFS:
        platform, run, sample = each.split("/")
        if platform == wc.platform and run == wc.run:
            for ref in SAMPLE_REFS[each]:
                ls.append(INTs + f"processed/cov_tmp/{platform}/{run}/{sample}__{ref}.finalCov.named.tsv")
    assert len(ls) >= 1
    return ls

rule cov__combine_cov_viz:
    input:
        mk_cov_input,
    output:
        cov = INTs + "cov/{platform}--{run}." + SITE_NAME + ".{ref}__all.finalCov.tsv",
    conda:
        REQS + "cat.conda.env.yaml"
    shell:
        "python-cath {input} {output}"


rule cov__plot_depths_for_report:
    "ngs core report"
    input:
        cov=INTs + "cov/{platform}--{run}." + SITE_NAME + ".{ref}__all.finalCov.tsv",
    output:
        o=SITE + "{platform}--{run}--{ref}__coverage.png",
    container:
        "docker://samesense/tidyverse-r:1"
    shell:
        "Rscript bin/plot_cov.R {wildcards.ref} {wildcards.run} {input} {output.o}"

def cov_plots_output(wc):
    ls = []
    for pr in RUN_REFS:
        platform, run = pr.split('--')
        refs = RUN_REFS[pr]
        for ref in refs:
            ls.append(SITE + f"{platform}--{run}--{ref}__coverage.png")
    assert len(ls) >= 1
    return ls

rule cov__target:
    input:
        cov_plots_output
