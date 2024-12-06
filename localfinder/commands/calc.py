# File: commands/calc.py

import sys
import json
import os
import pandas as pd
from localfinder.utils import (
    localPearson_and_enrichmentSignificance,
    localWeightedPearson_and_enrichmentSignificance,
    localSpearman_and_enrichmentSignificance,
    localWeightedSpearman_and_enrichmentSignificance,
    localMI_and_enrichmentSignificance,
    get_chromosomes_from_chrom_sizes
)

def main(args):
    track1_file = args.track1
    track2_file = args.track2
    output_dir = args.output_dir
    method = args.method
    method_params = args.method_params
    bin_number_of_window = args.bin_num
    step = args.step
    chroms = args.chroms
    chrom_sizes = args.chrom_sizes  # **Assuming chrom_sizes is now passed to calc.py**

    os.makedirs(output_dir, exist_ok=True)

    # **Modification Start**
    # If chroms is 'all', retrieve all chromosomes from chrom_sizes
    if chroms == ['all'] or chroms is None:
        chroms = get_chromosomes_from_chrom_sizes(chrom_sizes)
        print(f"'chroms' set to all chromosomes from chrom_sizes: {chroms}")
    else:
        print(f"'chroms' set to specified chromosomes: {chroms}")
    # **Modification End**

    # Read the binned tracks
    try:
        df1 = pd.read_csv(track1_file, sep='\t', header=None, names=['chr', 'start', 'end', 'readNum_1'])
        df2 = pd.read_csv(track2_file, sep='\t', header=None, names=['chr', 'start', 'end', 'readNum_2'])
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Merge the dataframes
    df = pd.merge(df1, df2, on=['chr', 'start', 'end'], how='inner')

    # # Remove bins with zero counts in both tracks
    # df = df[(df['readNum_1'] > 0) | (df['readNum_2'] > 0)]

    # Parse method_params if it's a string
    if isinstance(method_params, str):
        try:
            method_params = json.loads(method_params)
        except json.JSONDecodeError as e:
            print(f"Error parsing method_params JSON: {e}")
            sys.exit(1)

    # Call the selected method based on the alias
    if method == 'locP_and_ES':
        localPearson_and_enrichmentSignificance(
            df,
            bin_number_of_window=bin_number_of_window,
            step=step,
            output_dir=output_dir,
            chroms=chroms,
            chrom_sizes=chrom_sizes,
            **method_params
        )
    elif method == 'locWP_and_ES':
        localWeightedPearson_and_enrichmentSignificance(
            df,
            bin_number_of_window=bin_number_of_window,
            step=step,
            output_dir=output_dir,
            chroms=chroms,
            chrom_sizes=chrom_sizes,
            **method_params
        )
    elif method == 'locS_and_ES':
        localSpearman_and_enrichmentSignificance(
            df,
            bin_number_of_window=bin_number_of_window,
            step=step,
            output_dir=output_dir,
            chroms=chroms,
            chrom_sizes=chrom_sizes,
            **method_params
        )
    elif method == 'locWS_and_ES':
        localWeightedSpearman_and_enrichmentSignificance(
            df,
            bin_number_of_window=bin_number_of_window,
            step=step,
            output_dir=output_dir,
            chroms=chroms,
            chrom_sizes=chrom_sizes,
            **method_params
        )
    elif method == 'locMI_and_ES':
        localMI_and_enrichmentSignificance(
            df,
            bin_number_of_window=bin_number_of_window,
            step=step,
            output_dir=output_dir,
            chroms=chroms,
            chrom_sizes=chrom_sizes,
            **method_params
        )
    else:
        print(f"Unsupported method: {method}")
        sys.exit(1)
