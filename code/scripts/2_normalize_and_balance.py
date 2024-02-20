from typing import Any, TYPE_CHECKING
if TYPE_CHECKING:
    snakemake: Any = None

import pathlib
import subprocess
import cooler

resolutions = snakemake.config['resolutions']
labels = list(snakemake.config['labels'])

output_dir = snakemake.params['output_dir']
pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

input_dir = snakemake.params['input_dir']
hic_pro_results_dir = f'{input_dir}/hic-pro/hic_results/matrix/'

for resolution in resolutions:
    print(f'Normalizing and balancing at {resolution}bp resolution')
    for label in labels:
        bins_path = f'{hic_pro_results_dir}/{label}/raw/{resolution}/{label}_{resolution}_abs.bed'
        pixels_path = f'{hic_pro_results_dir}/{label}/raw/{resolution}/{label}_{resolution}.matrix'
        cool_path = f'{output_dir}/{label}_{resolution}_raw.cool'
        subprocess.run([
            'cooler', 'load',
            '--format', 'coo',
            '--assembly', 'mm10',
            bins_path,
            pixels_path,
            cool_path
        ], check=True)

    in_coolers = []
    out_coolers = []
    for label in labels:
        in_coolers.append(f'{output_dir}/{label}_{resolution}_raw.cool')
        out_coolers.append(f'{output_dir}/{label}_{resolution}_normalized.cool')
    # in_coolers = ' '.join(in_coolers)
    # out_coolers = ' '.join(out_coolers)
    subprocess.run([
        'hicNormalize',
        '-m', *in_coolers,
        '--normalize', 'smallest',
        '-o', *out_coolers
    ], check=True)

    for label in labels:
        in_cooler = f'{output_dir}/{label}_{resolution}_normalized.cool'
        out_cooler = f'{output_dir}/{label}_{resolution}_balanced.cool'
        subprocess.run([
            'cp', in_cooler, out_cooler
        ], check=True)
        cool = cooler.Cooler(out_cooler)
        balanced = cooler.balance_cooler(cool, store=True, cis_only=False, trans_only=False)
