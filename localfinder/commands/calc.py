# File: commands/calc.py

import os, sys, pandas as pd
from localfinder.utils import (
    locCor_and_ES,
    get_chromosomes_from_chrom_sizes
)

def main(args):
    track1_file = args.track1
    track2_file = args.track2
    output_dir = args.output_dir
    method = args.method
    percentile = args.percentile
    bin_number_of_window = args.binNum_window
    bin_number_of_peak = args.binNum_peak
    FC_thresh = args.FC_thresh
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


    # decide which correlation to use
    if method == 'locP_and_ES':
        corr_method = 'pearson'
    elif method == 'locS_and_ES':
        corr_method = 'spearman'
    else:
        print(f"Unsupported method: {method}")
        sys.exit(1)

        # ── RUN the analysis ──────────────────────────────────────────────────────
    locCor_and_ES(
        df,
        bin_number_of_window = bin_number_of_window,
        step                 = step,
        percentile           = percentile,
        FC_thresh            = FC_thresh,
        bin_number_of_peak   = bin_number_of_peak,
        corr_method          = corr_method,
        output_dir           = output_dir,
        chroms               = chroms
    )