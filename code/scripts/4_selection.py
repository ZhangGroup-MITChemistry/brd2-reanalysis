from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import subprocess
import bioframe as bf
import numpy as np
import scipy.stats

resolution = int(snakemake.wildcards['resolution'])
labels = list(snakemake.config['labels'])

eigen_output_dir = f"{snakemake.params['eigen_output_dir']}/{resolution}"

chromosomes = [f'chr{i}' for i in range(1, 20)]
df = bf.binnify(bf.read_chromsizes(f'{snakemake.config["MM10_DATA"]}/mm10.chrom.sizes'), resolution)

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

# use mm10.GCpt.tss.bedGraph to manually select and orient from both tsvs
phasing_tracks = bf.read_table(f'{MM10_DATA_CP}/mm10.GCpt.tss.bedGraph', names=['chrom', 'start', 'end', 'GCpt', 'tss'])
eigenvectors_cooltools = bf.read_table(f'{eigen_output_dir}/eigenvectors_cooltools.tsv', header=0)
eigenvectors_dchic = bf.read_table(f'{eigen_output_dir}/eigenvectors_dchic.tsv', header=0)

for method in ['cooltools', 'dchic']:
    egvec_table = eval(f'eigenvectors_{method}')
    compartments_table = df.copy()
    compartments_table[[f'{label}_compartments' for label in labels]] = np.nan
    for chrom in chromosomes:
        gc_track = phasing_tracks[phasing_tracks['chrom'] == chrom]['GCpt'].to_numpy()
        tss_track = phasing_tracks[phasing_tracks['chrom'] == chrom]['tss'].to_numpy()
        for label in labels:
            E1 = egvec_table[egvec_table['chrom'] == chrom][f'{label}_E1'].to_numpy()
            E2 = egvec_table[egvec_table['chrom'] == chrom][f'{label}_E2'].to_numpy()
            gc_corr_E1 = scipy.stats.spearmanr(gc_track, E1, nan_policy='omit').statistic
            gc_corr_E2 = scipy.stats.spearmanr(gc_track, E2, nan_policy='omit').statistic
            tss_corr_E1 = scipy.stats.spearmanr(tss_track, E1, nan_policy='omit').statistic
            tss_corr_E2 = scipy.stats.spearmanr(tss_track, E2, nan_policy='omit').statistic
            abs_corr_E1 = abs(gc_corr_E1) + abs(tss_corr_E1)
            abs_corr_E2 = abs(gc_corr_E2) + abs(tss_corr_E2)
            if abs_corr_E1 > abs_corr_E2:
                compartments_table.loc[compartments_table['chrom'] == chrom, f'{label}_compartments'] = np.sign(gc_corr_E1) * E1
            else:
                compartments_table.loc[compartments_table['chrom'] == chrom, f'{label}_compartments'] = np.sign(gc_corr_E2) * E2
    compartments_table.to_csv(snakemake.output[f'eigenvectors_{method}'], sep='\t', index=False, header=True)
