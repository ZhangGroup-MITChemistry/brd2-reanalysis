#!/bin/bash

#SBATCH -c 1
#SBATCH -o /home/gridsan/aathreya/files/projects/brd2/main/results/logs/%j.out

module purge
eval "$(conda shell.bash hook)"
conda activate dchic

snakemake --profile=profile
