rule deploy_coverage_app_multiSample:
    input:
        INTs + "flags/rsconnect",
        cov = INTs + "cov/{run}." + SITE_NAME + ".{ref}__all.finalCov.tsv",
        regions = INTs + "regions.standardized.tsv",
    output:
        ip_file = INTs + "dash-cov-ms/ip_files/{run}--{ref}--cov.ip",
    params:
        reqs = REQS + "dash-cov-requirements.txt",
        app_json = CONFIG + "dash-cov-json/{run}/{ref}/" + SITE_NAME,
        appdir = INTs + "dash-cov/tmp/app/{run}/{ref}/" + SITE_NAME + "/",
        app_template = "./bin/app_template_samples.py",
    conda:
        REQS + "rsconnect.conda.env.yaml"
    threads: 1
    shell:
        '''
        python bin/deploy_dash_cov_app_ms.py {wildcards.ref} {wildcards.run} \
        {input.cov} {input.regions} {params.app_template} {params.reqs} \
        {params.app_json} {params.appdir} {output.ip_file} "{wildcards.run}.{SITE_NAME}"
        '''


def mk_cov_vis_ls_for_run(wc):
    ls = []
    for ref in RUN_REFS[f"shortRead--{wc.run}"]:
        af = INTs + f"dash-cov-ms/ip_files/shortRead--{wc.run}--{ref}--cov.ip"
        ls.append(af)
    return ls
