from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import pathlib
import subprocess

root = snakemake.config['root']

RESOLUTIONS = list(snakemake.config['resolutions'])
RESOLUTIONS_STR = ' '.join([str(int(x)) for x in RESOLUTIONS])

config_hicpro = f"""# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = 48
SORT_RAM = 180000M
LOGFILE = hicpro.log

JOB_NAME = hic-pro
JOB_MEM = 4000M
JOB_WALLTIME = 24:00:00
JOB_QUEUE = xeon-p8
JOB_MAIL = ""

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _1
PAIR2_EXT = _2

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10

BOWTIE2_IDX_PATH = {root}/pkgs/HiC-Pro_3.1.0/annotation/Mouse/mm10
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive-local --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive-local --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = mm10
GENOME_SIZE = {root}/pkgs/HiC-Pro_3.1.0/annotation/chrom_mm10.sizes

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Capture Hi-C analysis
#######################################################################

CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = 
LIGATION_SITE = 
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = {RESOLUTIONS_STR}
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
"""

output_dir = snakemake.output['output_directory']
pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
with open(f'{output_dir}/config_hicpro.txt', 'w') as f:
    f.write(config_hicpro)

cmd0 = f'rm -rf {output_dir}/hic-pro'
print(cmd0)
shout = subprocess.Popen(cmd0, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)

cmd1 = f'{root}/pkgs/HiC-Pro_3.1.0/bin/HiC-Pro -i {snakemake.config["SRA_DATA_DIR"]} -o {output_dir}/hic-pro -c {output_dir}/config_hicpro.txt -p'
print(cmd1)
shout = subprocess.Popen(cmd1, shell=True, cwd=output_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)

cmd2 = f'sbatch HiCPro_step1_hic-pro.sh'
print(cmd2)
shout = subprocess.Popen(cmd2, shell=True, cwd=f'{output_dir}/hic-pro', stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)
step1_jobid = shout.split()[-1]

cmd3 = f'sbatch --dependency=afterok:{step1_jobid} --wait HiCPro_step2_hic-pro.sh'
print(cmd3)
shout = subprocess.Popen(cmd3, shell=True, cwd=f'{output_dir}/hic-pro', stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)
