#!/usr/bin/env python
import sys, os, uuid, requests


def stage_app_file(tmp_dir, app_template, ref, regionTsv, sitename):
    regionFile = regionTsv.split("/")[-1]
    with open(tmp_dir + "app.py", "w") as fout:
        with open(app_template) as f:
            for line in f:
                if "{replace_ref" in line:
                    ll = line.rstrip().replace("{replace_ref}", ref)
                elif "{replace_site" in line:
                    ll = line.rstrip().replace("{replace_sitename}", sitename)
                elif "{regionsStandardizedTsv}" in line:
                    ll = line.rstrip().replace("{regionsStandardizedTsv}", regionFile)
                else:
                    ll = line.rstrip()
                print(ll, file=fout)


def stage(
    tmp_dir, cov_tsv, region_tsv, app_json, app_template, ref, u, dash_req_py, sitename
):
    os.system(f"mkdir -p {tmp_dir}")
    os.system(f"cp {cov_tsv} {tmp_dir}")
    os.system(f"cp {region_tsv} {tmp_dir}")
    os.system(f"cp {dash_req_py} {tmp_dir}requirements.txt")
    stage_app_file(tmp_dir, app_template, ref, region_tsv, sitename)
    if os.path.exists(app_json):
        d = f"{tmp_dir}rsconnect-python/"
        os.system(f"mkdir -p {d}")
        os.system(f"cp {app_json} {d}{u}.json")


def set_access_all(url, domain, api_endpoint, headers):
    """Set access to 'all' for content published to a url
    Returns guid for that content
    """
    access_all = {"access_type": "all"}
    guid = url.split("/")[-2]
    api_url = "".join([domain, api_endpoint, "experimental/content/", guid])
    response = requests.post(api_url, headers=headers, json=access_all)
    return guid


def deploy(tmp_dir, ip_log):
    cmd = f"cd {tmp_dir}; rsconnect deploy dash . > {ip_log}.tmp "
    os.system(cmd)
    cmd = f"grep Direct {ip_log}.tmp | grep URL | cut -f 8 -d ' ' > {ip_log}"
    os.system(cmd)

    with open(ip_log) as f:
        cov_url = f.readline().strip()
    DOMAIN = "https://connect.sparkds.io"
    API_ENDPOINT = "/__api__/v1/"
    API_KEY = "za4RPQdH58QuGHPTZbHvU087nF4TlkD3"
    HEADER = {"Authorization": "Key " + API_KEY}

    set_access_all(cov_url, DOMAIN, API_ENDPOINT, HEADER)


def main():
    (
        ref,
        minusBamStr,
        cov_tsv,
        region_tsv,
        app_template_py,
        dash_req_py,
        app_json,
        app_dir,
        ip_outfile,
        sitename,
    ) = sys.argv[1:]
    u = str(uuid.uuid3(uuid.NAMESPACE_URL, f"{ref}__{minusBamStr}-cov")) + "-cov"
    tmp_dir = f"{app_dir}{u}/"
    # make app dir and file
    # use app json if available
    stage(
        tmp_dir,
        cov_tsv,
        region_tsv,
        app_json,
        app_template_py,
        ref,
        u,
        dash_req_py,
        sitename,
    )

    # deploy
    deploy(tmp_dir, ip_outfile)

    # cp dash json for next time
    # os.system(f"cp {tmp_dir}rsconnect-python/{u}.json {app_json}")

    # delete tmp dir
    # os.system(f'rm -rf {tmp_dir}')


if __name__ == "__main__":
    main()
