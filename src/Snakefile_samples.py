from dataclasses import dataclass
from collections import defaultdict
import glob
import gzip
import csv
import re

@dataclass
class seqsample:
    '''
    A11-15_S12_L001, miseq
    s3://s3-vendor-device-rd-1250-mseq/211209_M05240_0015_000000000-JVBHH/Alignment_1/20211211_152909/Fastq/A11-15_S12_L001_R1_001.fastq.gz
    A11-15_S12_L001_R1_001
    run: 211209_M05240_0015_000000000-JVBHH
    '''
    name: str
    r1_path: str
    r1_name: str
    r2_path: str
    r2_name: str
    platform: str
    run: str
    alignment: str


def fq_has_reads(fq):
    reads = 0
    with gzip.open(fq, "rt") as handle:
        for line in handle:
            reads += 1
            if reads == 40:
                break
    return reads > 3

def parse_seq_bucket():
    '''sample names: run/platform/name'''
    platform = 'shortRead'
    ret_samples = {}
    sample_reads = defaultdict(dict)

    data_dir = DATA + 'raw/shortRead/'
    globs = glob.glob(data_dir + '**/Alignment_*/**/Fastq/*.fastq.gz') + glob.glob(data_dir + '**/Analysis/1/Data/fastq/*.fastq.gz')
        
    for fq in globs:
        s = str(fq)
        # miseq data structure Alignment_1/**/Fastq/*.fastq.gz
        if re.search('Alignment', s):
            sample = s.split('_L001')[0].split('/')[-1]
            run = s.split('/')[-5]
            a = s.split('/Alignment_')[-1].split('/')[0]
            if fq_has_reads(fq):
                sample = f'{platform}/{run}/{sample}_a{a}'
                # A11-15_S12_L001_R1_001.fastq.gz
                r = s.split('_L001_')[1].split('_')[0]
                assert r in ('R1', 'R2')
                sample_reads[sample][r] = (str(fq), s.split('/')[-1].split('.fastq')[0])
        # next-seq data structure Analysis/1/Data/fastq/*.fastq.gz
        elif re.search('Analysis', s):
            if '_001' in s:
                sample = s.split('_001')[0].split('/')[-1][:-3]
            else:
                sample = s.split('.')[0].split('/')[-1][:-3]
            # print(sample)
            run = s.split('/')[-6]
            a = '1' #s.split('/Alignment_')[-1].split('/')[0]
            if fq_has_reads(fq):
                sample = f'{platform}/{run}/{sample}_a{a}'
                if '_001' in s:
                    r = s.split('_001')[0].split('_')[-1]
                else:
                    r = s.split('.')[0].split('_')[-1]
                assert r in ('R1', 'R2')
                sample_reads[sample][r] = (str(fq), s.split('/')[-1].split('.fastq')[0])


    for sample in sample_reads:
        #print(sample_reads[sample])
        r1_path, r1_name = sample_reads[sample]['R1']
        if 'R2' in sample_reads[sample]:
            r2_path, r2_name = sample_reads[sample]['R2']
        else:
            r2_path, r2_name = "", ""

        if 'L001' in r1_path:
            run = r1_path.split('/')[-5]
        else:
            run = r1_path.split('/')[-6]

        if all(keyword not in sample for keyword in ['deter', 'unknown']):
            ret_samples[sample] = seqsample(sample, r1_path, r1_name, r2_path, r2_name, 'shortRead', run, sample[-1])
    return ret_samples

def find_samples():
    return parse_seq_bucket()

def get_bucket_samples():
    samples = {}
    sd = find_samples()
    for sample in sd:
        samples[sample] = sd[sample]
    return samples

def list_samples():
    sample_obj_dict = get_bucket_samples()
    samples = []
    for sample in sample_obj_dict:
        samples.append(sample)
    return samples, sample_obj_dict

def get_platform_runs():
    platform_runs = defaultdict(dict)
    ls = []
    for sample in SAMPLES_PRE:
        platform, run, name = sample.split("/")
        platform_runs[f"{platform}--{run}"][sample] = True
    return platform_runs

SAMPLES_PRE, SAMPLE_OBJ_DICT = list_samples()
# print(SAMPLES_PRE)
SAMPLES = dict()
RUNS_PRE = get_platform_runs()
# print(RUNS_PRE)
RUNS = {}
SAMPLE_REFS = defaultdict(dict)
RUN_REFS = defaultdict(dict)
REFS = {}
RUN_HOST = {}
RUN_LIB_TYPE = {}
RUN_REP = defaultdict(list)
with open(MAP_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        RUN_HOST[f'{row["run"]}'] = f'{row["host_ref"]}'
        RUN_LIB_TYPE[f'{row["run"]}'] = f'{row["library_type"]}'
        rep = f'{row["tech_replicate"]}'
        if rep not in RUN_REP[f'{row["run"]}']:
            RUN_REP[f'{row["run"]}'].append(rep)
        key = f'shortRead/{row["run"]}/{row["Sample"]}'
        for ref in row["reference"].split(";"):
            if True:  # not "210907" in key:
                SAMPLES[key] = True  # SAMPLES_PRE[key]
                SAMPLE_REFS[key][ref] = True
                run_key = f'shortRead--{row["run"]}'
                RUN_REFS[run_key][ref] = True
                RUNS[run_key] = RUNS_PRE[run_key]
                REFS[ref] = True
                
# print(RUNS)
# print(REFS)
# print(RUN_HOST)
# print(RUN_LIB_TYPE)
# print(RUN_REP)
# print(RUNS_PRE)
# print(SAMPLES)
# print(SAMPLES_PRE)
# print(RUN_REFS) 
# print(SAMPLE_REFS)
# print(SAMPLE_OBJ_DICT)
