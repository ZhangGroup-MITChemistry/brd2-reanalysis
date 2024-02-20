from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import subprocess
# import pathlib
# import cooler
# import cooltools
# import bioframe as bf
# import pandas as pd
# import numpy as np

resolution = int(snakemake.wildcards['resolution'])
labels = list(snakemake.config['labels'])

eigen_output_dir = f"{snakemake.params['eigen_output_dir']}/{resolution}"

MM10_DATA_CP = f'{snakemake.config["MM10_DATA"]}_{resolution}_goldenpathData'
cmd0 = f'cp -r {snakemake.config["MM10_DATA"]} {MM10_DATA_CP}'
print(cmd0)
shout = subprocess.Popen(cmd0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)

# run dcHiC select
cmd1 = f'conda run -n dchic Rscript {snakemake.config["DCHIC_DIR"]}/dchicf.r --file dchic_input.txt --pcatype select --dirovwt T --genome mm10 --gfolder {MM10_DATA_CP}'
print(cmd1)
shout = subprocess.Popen(cmd1, shell=True, cwd=f'{eigen_output_dir}/dchic/analysis', stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)

# use GCpt.tss.bedgraph to manually select and orient from both tsvs

# return a bedgraph of the compartments genomewide

for _, file in snakemake.output.items():
    print(file)
    with open(file, 'w') as f:
        f.write('')
