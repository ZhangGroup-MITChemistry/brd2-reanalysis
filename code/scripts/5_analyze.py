from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import bioframe as bf
import numpy as np
import cooler
import cooltools

resolution = int(snakemake.wildcards['resolution'])
labels = list(snakemake.config['labels'])
chromosomes = [f'chr{i}' for i in range(1, 20)]

cool_dir = f"{snakemake.params['cool_input_dir']}"

# saddle plots
q = 1
Q_LO, Q_HI = q / 100, 1 - q / 100
N_GROUPS = int(100 / q - 2)

MM10_DATA_CP = f'{snakemake.config["MM10_DATA"]}_{resolution}_goldenpathData'

# reference_track = bf.read_table(f'{MM10_DATA_CP}/mm10.GCpt.bedGraph', names=['chrom', 'start', 'end', 'GCpt'])

for method in ['cooltools', 'dchic']:
    egvec_table = bf.read_table(snakemake.input[f'eigenvectors_{method}'], header=0)
    reference_track = egvec_table[['chrom', 'start', 'end', 'Control_B2_compartments']].copy()

    results = {}

    for label in labels:
        results[label] = {'whole_genome': {}, 'per_chromosome': {}}
        
        cool = cooler.Cooler(f'{cool_dir}/{label}_{resolution}_balanced.cool')
        expected = cooltools.expected_cis(cool)
        interaction_sum, interaction_count =  cooltools.saddle(
                cool,
                expected,
                reference_track,
                'cis',
                n_bins=N_GROUPS,
                qrange=(Q_LO,Q_HI),
        )

        digitized_track, _ = cooltools.digitize(reference_track, N_GROUPS, qrange=(Q_LO,Q_HI))

        C = (interaction_sum / interaction_count)[1:-1, 1:-1]
        groupmean = reference_track[reference_track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()[1:-1]
        saddle_strength = cooltools.api.saddle.saddle_strength(interaction_sum, interaction_count)[1:]

        results[label]['whole_genome']['saddle_grid'] = C.copy()
        results[label]['whole_genome']['saddle_hist'] = groupmean.copy()
        results[label]['whole_genome']['saddle_strength'] = saddle_strength.copy()

        for chrom in chromosomes:
            reference_track_chrom = reference_track[reference_track['chrom'] == chrom]
            interaction_sum, interaction_count =  cooltools.saddle(
                    cool,
                    expected,
                    reference_track_chrom,
                    'cis',
                    n_bins=N_GROUPS,
                    qrange=(Q_LO,Q_HI),
            )

            digitized_track, _ = cooltools.digitize(reference_track_chrom, N_GROUPS, qrange=(Q_LO,Q_HI))

            C = (interaction_sum / interaction_count)[1:-1, 1:-1]
            groupmean = reference_track_chrom[reference_track_chrom.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()[1:-1]
            saddle_strength = cooltools.api.saddle.saddle_strength(interaction_sum, interaction_count)[1:]
            
            results[label]['per_chromosome'][chrom] = {'saddle_grid': C.copy(), 'saddle_hist': groupmean.copy(), 'saddle_strength': saddle_strength.copy()}

    np.savez(snakemake.output[f'saddles_{method}'], **results)
