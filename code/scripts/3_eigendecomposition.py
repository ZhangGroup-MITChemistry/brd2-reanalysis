from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import pathlib
import subprocess
import cooler
import cooltools
import bioframe as bf
import pandas as pd
import numpy as np

resolution = int(snakemake.wildcards['resolution'])
labels = list(snakemake.config['labels'])

input_dir = snakemake.params['input_dir']
output_dir = f"{snakemake.params['output_dir']}/{resolution}"
pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

chromosomes = [f'chr{i}' for i in range(1, 20)]
df = bf.binnify(bf.read_chromsizes(f'{snakemake.config["MM10_DATA"]}/mm10.chrom.sizes'), resolution)

# dcHiC
pathlib.Path(f'{output_dir}/dchic/input').mkdir(parents=True, exist_ok=True)
for label in labels:
    cool = cooler.Cooler(f'{input_dir}/{label}_{resolution}_balanced.cool')
    cool.matrix(balance=True, as_pixels=True)[:][['bin1_id', 'bin2_id', 'balanced']].to_csv(f'{output_dir}/dchic/input/{label}_{resolution}_balanced.matrix', sep='\t', index=False, header=False)
    bins = cool.bins()[:].reset_index()[['chrom', 'start', 'end', 'index']]
    bins = bins[bins['chrom'].isin(chromosomes)]
    bins.to_csv(f'{output_dir}/dchic/input/{label}_{resolution}_balanced_bins.bed', sep='\t', index=False, header=False)

# FIXME
samples = {
    'Control_B2': 'Control_B2',
    'Brd2_dep_rep1': 'Brd2_dep',
    'Brd2_dep_rep2': 'Brd2_dep',
}
# FIXME END

pathlib.Path(f'{output_dir}/dchic/analysis').mkdir(parents=True, exist_ok=True)
with open(f'{output_dir}/dchic/analysis/dchic_input.txt', 'w') as f:
    s = ''
    for label in labels:
        BINS_PATH = f'{output_dir}/dchic/input/{label}_{resolution}_balanced_bins.bed'
        PIXELS_PATH = f'{output_dir}/dchic/input/{label}_{resolution}_balanced.matrix'
        s += f'{PIXELS_PATH}\t{BINS_PATH}\t{label}\t{samples[label]}\n'
    f.write(s)
  
cmd1 = f'conda run -n dchic Rscript {snakemake.config["DCHIC_DIR"]}/dchicf.r --file dchic_input.txt --pcatype cis --dirovwt T &> {output_dir}/dchic/analysis/dchicf_cis.log'
print(cmd1)
shout = subprocess.Popen(cmd1, shell=True, cwd=f'{output_dir}/dchic/analysis', stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read().decode('utf-8')
print(shout)

eigenvectors_dchic_df = df.copy()
for label in labels:
    all_chrom_PCs = []
    for chrom in chromosomes:
        chrom_PCs = bf.read_table(f'{output_dir}/dchic/analysis/{label}_pca/intra_pca/{label}_mat/{chrom}.pc.txt', header=0).rename(columns={'chr': 'chrom', 'PC1': f'{label}_E1', 'PC2': f'{label}_E2'}).drop(columns=['index'])
        all_chrom_PCs.append(chrom_PCs.copy())
    all_chrom_PCs = pd.concat(all_chrom_PCs)
    eigenvectors_dchic_df = eigenvectors_dchic_df.merge(all_chrom_PCs, on=['chrom', 'start', 'end'], how='left')

eigenvectors_dchic_df.to_csv(snakemake.output['eigenvectors_dchic'], sep='\t', index=False)

# cooltools
eigenvectors_cooltools_df = df.copy()
for label in labels:
    cool = cooler.Cooler(f'{input_dir}/{label}_{resolution}_balanced.cool')
    eigs = cooltools.eigs_cis(cool, n_eigs=2)
    eigenvectors = eigs[1][['chrom', 'start', 'end', 'E1', 'E2']].rename(columns={'E1': f'{label}_E1', 'E2': f'{label}_E2'})
    eigenvectors_cooltools_df = eigenvectors_cooltools_df.merge(eigenvectors, on=['chrom', 'start', 'end'], how='left')

eigenvectors_cooltools_df.to_csv(snakemake.output['eigenvectors_cooltools'], sep='\t', index=False)
