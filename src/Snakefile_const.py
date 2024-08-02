import pathlib
D = config['PROJ_DIR_v100'] 
DATA = D + 'data/'
INTs = DATA + 'interim/'
REQS = D + 'reqs/'
CONFIG = D + 'configs/'

TRIM_FQ_DIR = config['TRIM_FQ_DIR']
SITE = config['SITE']
SITE_NAME = config['SITE_NAME']

# MAP_FILE = CONFIG + "master_mapping.tsv"
MAP_FILE = CONFIG + "master_mapping_anno.tsv"
# MAP_FILE = CONFIG + "master_mapping_rerun.tsv"
QC = config['QC']
MERGE_TAB = config['MERGE_TAB']
HOST_REF = ['hg38', 'mm39']
# print(INTs)