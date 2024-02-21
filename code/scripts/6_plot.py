from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import pathlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm

resolution = int(snakemake.wildcards['resolution'])
labels = list(snakemake.config['labels'])
chromosomes = [f'chr{i}' for i in range(1, 20)]

output_dir = snakemake.output['plots']

compartments_cooltools_df = pd.read_csv(snakemake.input['compartments_cooltools'], sep='\t', header=0)
compartments_dchic_df = pd.read_csv(snakemake.input['compartments_dchic'], sep='\t', header=0)

compartments_cooltools_df = compartments_cooltools_df.rename(columns={f'{x}': f'{x}_cooltools' for x in compartments_cooltools_df.columns if x not in ['chrom', 'start', 'end']})
compartments_dchic_df = compartments_dchic_df.rename(columns={f'{x}': f'{x}_dchic' for x in compartments_dchic_df.columns if x not in ['chrom', 'start', 'end']})

results_df = compartments_cooltools_df.merge(compartments_dchic_df, on=['chrom', 'start', 'end'], how='left')

control_sample = 'Control_B2'
expt = 'Brd2_dep'
replicates = [f'{expt}_rep1', f'{expt}_rep2']

plots = []

for method in ['cooltools', 'dchic']:
    for rep in range(len(replicates)):
        plots.append(
            {
                'Control': f'{control_sample}_compartments_{method}',
                'Brd2 dep': f'{replicates[rep]}_compartments_{method}',
                'method': method,
                'replicate': f'Replicate {rep + 1}',
            }
        )

output_folders = {
    'whole_genome': ['kdeplots', 'scatterplots', 'saddleplots', 'saddle_strengths'],
    'per_chromosome': ['kdeplots', 'scatterplots', 'profiles', 'saddleplots', 'saddle_strengths'],
}

for f1, f2 in output_folders.items():
    for f in f2:
        pathlib.Path(f'{output_dir}/{f1}/{f}').mkdir(parents=True, exist_ok=True)

for plot in plots:
    # plot kdeplot of control vs expt for whole genome
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
    sns.kdeplot(data=results_df[results_df['chrom'].isin(chromosomes)], x=plot['Control'], alpha=0.25, label='Control', ax=ax, color='k', fill=True, linewidth=1)
    sns.kdeplot(data=results_df[results_df['chrom'].isin(chromosomes)], x=plot['Brd2 dep'], alpha=0.25, label='Brd2 dep', ax=ax, color='xkcd:mint green', fill=True, linewidth=1, edgecolor='k')
    if plot['method'] == 'cooltools':
        ax.set_xlim(-2, 2)
    ax.set(xlabel='Compartment scores', ylabel='Density', title=f'{plot["method"]} {plot["replicate"]} {resolution / 1000:0.0f}kb')
    ax.legend(loc='upper right')
    fig.tight_layout()
    plt.savefig(f'{output_dir}/whole_genome/kdeplots/{plot["method"]}_{plot["replicate"]}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot kdeplot of control vs expt for each chromosome
    fig, axes = plt.subplots(4, 5, figsize=(24, 16), dpi=300, sharex=True, sharey=True)
    for ax in axes.flatten():
        ax.set_visible(False)
    for i, chrom in enumerate(chromosomes):
        select_df = results_df[results_df['chrom']==chrom]
        ax = axes.flatten()[i]
        ax.set_visible(True)
        sns.kdeplot(data=select_df, x=plot['Control'], alpha=0.25, label='Control', ax=ax, color='k', fill=True, linewidth=1)
        sns.kdeplot(data=select_df, x=plot['Brd2 dep'], alpha=0.25, label='Brd2 dep', ax=ax, color='xkcd:mint green', fill=True, linewidth=1, edgecolor='k')
        if plot['method'] == 'cooltools':
            ax.set_xlim(-2, 2)
        ax.set(xlabel=None, ylabel=None, title=chrom)
    fig.legend(*ax.get_legend_handles_labels(), loc='center', bbox_to_anchor=(0.9, 0.125))
    fig.supxlabel('Compartment scores')
    fig.supylabel('Density')
    fig.suptitle(f'{plot["method"]} {plot["replicate"]} {resolution / 1000:0.0f}kb', fontsize=16)
    fig.tight_layout()
    plt.savefig(f'{output_dir}/per_chromosome/kdeplots/{plot["method"]}_{plot["replicate"]}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot scatter of control vs expt for whole genome
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
    x, y = results_df[plot['Control']], results_df[plot['Brd2 dep']]
    mask = ~x.isna() & ~y.isna()
    sd = np.sqrt(np.cov(x[mask], y[mask])[0, 1])

    ax.scatter(x, y, alpha=0.25, color='C2', s=1)
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    p0, p1 = min(x0, y0), max(x1, y1)
    ax.plot([p0, p1], [p0, p1], 'r')
    ax.plot([p0, p1], [p0 + sd, p1 + sd], 'b--')
    ax.plot([p0, p1], [p0 - sd, p1 - sd], 'b--')
    ax.set(xlabel='Control', ylabel='Brd2 dep', title=f'{plot["method"]} {plot["replicate"]} {resolution / 1000:0.0f}kb')
    fig.tight_layout()
    plt.savefig(f'{output_dir}/whole_genome/scatterplots/{plot["method"]}_{plot["replicate"]}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot scatter of control vs expt for each chromosome
    fig, axes = plt.subplots(4, 5, figsize=(24, 16), dpi=300, sharex=True, sharey=True)
    for ax in axes.flatten():
        ax.set_visible(False)
    for i, chrom in enumerate(chromosomes):
        select_df = results_df[results_df['chrom']==chrom]
        ax = axes.flatten()[i]
        ax.set_visible(True)
        
        x, y = select_df[plot['Control']], select_df[plot['Brd2 dep']]
        mask = ~x.isna() & ~y.isna()
        sd = np.sqrt(np.cov(x[mask], y[mask])[0, 1])

        ax.scatter(x, y, alpha=0.25, color='C2', s=1)
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        p0, p1 = min(x0, y0), max(x1, y1)
        ax.plot([p0, p1], [p0, p1], 'r')
        ax.plot([p0, p1], [p0 + sd, p1 + sd], 'b--')
        ax.plot([p0, p1], [p0 - sd, p1 - sd], 'b--')
        ax.set(xlabel=None, ylabel=None, title=chrom)
    fig.supxlabel('Control')
    fig.supylabel('Brd2 dep')
    fig.suptitle(f'{plot["method"]} {plot["replicate"]} {resolution / 1000:0.0f}kb', fontsize=16)
    fig.tight_layout()
    plt.savefig(f'{output_dir}/per_chromosome/scatterplots/{plot["method"]}_{plot["replicate"]}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot compartment profiles for each chromosome
    for i, chrom in enumerate(chromosomes):
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        select_df = results_df[results_df['chrom']==chrom]
        sns.lineplot(data=select_df, x='start', y=plot['Control'], label='Control', color='black', ax=ax)
        sns.lineplot(data=select_df, x='start', y=plot['Brd2 dep'], label='Brd2 dep', color='green', ax=ax)
        ax.axhline(0, color='black', lw=1)

        ax.set_xlabel(chrom)
        ax.set_ylabel('Compartment score')
        ax.set_title(f'{plot["method"]} {plot["replicate"]} {resolution / 1000:0.0f}kb', fontsize=16)
        fig.tight_layout()
        plt.savefig(f'{output_dir}/per_chromosome/profiles/{plot["method"]}_{plot["replicate"]}_{resolution / 1000:0.0f}kb_{chrom}.pdf', bbox_inches='tight', dpi=300)


q = 1
Q_LO, Q_HI = q / 100, 1 - q / 100
N_GROUPS = int(100 / q - 2)

binedges = np.linspace(Q_LO, Q_HI, N_GROUPS + 1)
X, Y = np.meshgrid(binedges, binedges)

for method in ['cooltools', 'dchic']:
    saddles = np.load(snakemake.input[f'saddles_{method}'], allow_pickle=True)

    # plot saddle plots for whole genome
    fig = plt.figure(figsize=(16, 6.5), dpi=300)
    subfigs = fig.subfigures(1, 3)

    for i, label in enumerate(labels):
        C = saddles[label].item()['whole_genome']['saddle_grid']
        groupmean = saddles[label].item()['whole_genome']['saddle_hist']

        subfig = subfigs[i]
        axes = subfig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0)

        gs = axes[0]
        ax0 = subfig.add_subplot(gs)
        im = ax0.pcolormesh(X, Y, C, cmap='seismic', norm=LogNorm(vmin=0.5, vmax=2))
        ax0.set_ylim(ax0.get_ylim()[1], ax0.get_ylim()[0])
        ax0.set_xlim(Q_LO, Q_HI)

        gs = axes[1]
        ax1 = subfig.add_subplot(gs)
        ax1.bar(binedges, width=1/len(binedges), height=groupmean, align="edge", edgecolor='black')
        ax1.set_xlim(Q_LO, Q_HI)
        ax1.yaxis.set_visible(False)

        subfig.suptitle(f'{label}')

    cax = subfig.add_axes([ax0.get_position().x1+0.1, ax0.get_position().y0, 0.05, ax0.get_position().height])
    cbar = subfig.colorbar(im, cax=cax, ticks=[0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2], format='{x:0.1f}', label='Average O/E contact frequency')

    fig.suptitle(f'Whole genome saddle plots at {resolution / 1000:0.0f}kb', fontsize=16)
    plt.savefig(f'{output_dir}/whole_genome/saddleplots/{method}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot saddle plots for each chromosome
    for chrom in chromosomes:
        fig = plt.figure(figsize=(16, 6.5), dpi=300)
        subfigs = fig.subfigures(1, 3)
        for i, label in enumerate(labels):
            C = saddles[label].item()['per_chromosome'][chrom]['saddle_grid']
            groupmean = saddles[label].item()['per_chromosome'][chrom]['saddle_hist']

            subfig = subfigs[i]
            axes = subfig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0)

            gs = axes[0]
            ax0 = subfig.add_subplot(gs)
            im = ax0.pcolormesh(X, Y, C, cmap='seismic', norm=LogNorm(vmin=0.5, vmax=2))
            ax0.set_ylim(ax0.get_ylim()[1], ax0.get_ylim()[0])
            ax0.set_xlim(Q_LO, Q_HI)

            gs = axes[1]
            ax1 = subfig.add_subplot(gs)
            ax1.bar(binedges, width=1/len(binedges), height=groupmean, align="edge", edgecolor='black')
            ax1.set_xlim(Q_LO, Q_HI)
            ax1.yaxis.set_visible(False)

            subfig.suptitle(f'{label}')

        cax = subfig.add_axes([ax0.get_position().x1+0.1, ax0.get_position().y0, 0.05, ax0.get_position().height])
        cbar = subfig.colorbar(im, cax=cax, ticks=[0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2], format='{x:0.1f}', label='Average O/E contact frequency')

        fig.suptitle(f'{chrom} saddle plots at {resolution / 1000:0.0f}kb', fontsize=16)
        plt.savefig(f'{output_dir}/per_chromosome/saddleplots/{method}_{chrom}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot saddle strengths for whole genome
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
    for label in labels:
        css = saddles[label].item()['whole_genome']['saddle_strength']
        ax.step(binedges, css, where='pre', label=label)
    ax.set_xlabel('extent')
    ax.set_ylabel('(AA + BB) / (AB + BA)')
    ax.set_title(f'Whole genome saddle strength profile at {resolution / 1000:0.0f}kb')
    ax.axhline(1, c='grey', lw=1)
    ax.legend()
    fig.tight_layout()
    plt.savefig(f'{output_dir}/whole_genome/saddle_strengths/{method}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)

    # plot saddle strengths for each chromosome
    for chrom in chromosomes:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=300)
        for label in labels:
            css = saddles[label].item()['per_chromosome'][chrom]['saddle_strength']
            ax.step(binedges, css, where='pre', label=label)
        ax.set_xlabel('extent')
        ax.set_ylabel('(AA + BB) / (AB + BA)')
        ax.set_title(f'{chrom} saddle strength profile at {resolution / 1000:0.0f}kb')
        ax.axhline(1, c='grey', lw=1)
        ax.legend()
        fig.tight_layout()
        plt.savefig(f'{output_dir}/per_chromosome/saddle_strengths/{method}_{chrom}_{resolution / 1000:0.0f}kb.pdf', bbox_inches='tight', dpi=300)
