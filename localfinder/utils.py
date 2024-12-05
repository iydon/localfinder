# File: utils.py
import argparse
import numpy as np
import os
import pandas as pd
import shutil
import subprocess
import sys
import tempfile

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr, norm, spearmanr
from scipy.special import digamma
from sklearn.feature_selection import mutual_info_regression
from sklearn.neighbors import NearestNeighbors


def check_external_tools():
    required_tools = ['bedtools', 'bigWigToBedGraph', 'samtools']
    missing_tools = []
    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)

    if missing_tools:
        print("Error: The following required external tools are not installed or not found in PATH:")
        for tool in missing_tools:
            print(f" - {tool}")
        print("Please install them and ensure they are available in your system PATH.")
        sys.exit(1)


def get_chromosomes_from_chrom_sizes(chrom_sizes_file):
    """
    Retrieves all chromosome names from the chromosome sizes file.

    Parameters:
    - chrom_sizes_file (str): Path to the chromosome sizes file.

    Returns:
    - chromosomes (list): List of chromosome names.
    """
    try:
        df_chrom = pd.read_csv(chrom_sizes_file, sep='\t', header=None, names=['chr', 'size'])
        chromosomes = df_chrom['chr'].tolist()
        return chromosomes
    except Exception as e:
        print(f"Error reading chromosome sizes file {chrom_sizes_file}: {e}")
        sys.exit(1)
        
def process_and_bin_file(input_file, output_file, bin_size, chrom_sizes, chroms=None):
    """
    Processes an input file by detecting its format, converting it to BedGraph if necessary,
    and binning it into fixed-size bins. All intermediate files are handled in a temporary directory.
    """
    # Create a temporary directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_bedgraph = os.path.join(temp_dir, 'temp.bedgraph')

        # Detect the format based on the file extension
        input_format = detect_format(input_file)
        if input_format is None:
            print(f"Could not determine the format of the file: {input_file}")
            sys.exit(1)

        try:
            # Convert input to BedGraph if necessary
            if input_format == 'bigwig':
                if chroms:
                    # Process each chromosome separately
                    temp_bedgraph_files = []
                    for chrom in chroms:
                        temp_bedgraph_chrom = os.path.join(temp_dir, f'temp_{chrom}.bedgraph')
                        # Check the chromosome format by printing
                        print(f"Processing chromosome: {chrom}")
                        cmd = ['bigWigToBedGraph', '-chrom=' + chrom, input_file, temp_bedgraph_chrom]
                        subprocess.run(cmd, check=True)
                        temp_bedgraph_files.append(temp_bedgraph_chrom)

                    # Concatenate the per-chromosome bedgraph files into one
                    with open(temp_bedgraph, 'w') as outfile:
                        for fname in temp_bedgraph_files:
                            with open(fname) as infile:
                                outfile.write(infile.read())
                else:
                    # Process the entire file if no chroms specified
                    cmd = ['bigWigToBedGraph', input_file, temp_bedgraph]
                    subprocess.run(cmd, check=True)

            elif input_format == 'bedgraph':
                # Copy input to temporary BedGraph file
                subprocess.run(['cp', input_file, temp_bedgraph], check=True)

            else:
                print(f"Unsupported format: {input_format}")
                sys.exit(1)

            # Ensure the BedGraph file contains both chromosomes
            with open(temp_bedgraph, 'r') as temp_file:
                temp_data = temp_file.readlines()

            # Check if both chromosomes appear in the file
            found_chroms = set()
            for line in temp_data:
                chrom = line.split('\t')[0]
                found_chroms.add(chrom)

            print(f"Chromosomes found in {input_file}: {found_chroms}")

            # Bin the BedGraph data
            bin_bedgraph(temp_bedgraph, output_file, bin_size, chrom_sizes, chroms)

        finally:
            # Temporary directory and files are automatically cleaned up
            pass


def detect_format(filename):
    """
    Infers the file format based on the file extension.

    Parameters:
    - filename (str): The name or path of the file.

    Returns:
    - format (str): The inferred format ('bam', 'sam', 'bedgraph', 'bigwig'), or None if unknown.
    """
    extension = os.path.splitext(filename)[1].lower()
    if extension == '.bam':
        return 'bam'
    elif extension == '.sam':
        return 'sam'
    elif extension in ['.bedgraph', '.bdg']:
        return 'bedgraph'
    elif extension in ['.bigwig', '.bw']:
        return 'bigwig'
    else:
        return None


def bin_bedgraph(input_bedgraph, output_bedgraph, bin_size, chrom_sizes, chroms=None):
    """
    Bins the BedGraph data into fixed-size bins, replacing missing values with zero.
    Intermediate files are stored in a temporary directory that is deleted after processing.

    Parameters:
    - input_bedgraph (str): Path to the input BedGraph file.
    - output_bedgraph (str): Path to the output binned BedGraph file.
    - bin_size (int): Size of the bins in base pairs.
    - chrom_sizes (str): Path to the chromosome sizes file.
    - chroms (list): List of chromosomes to process. If None, all chromosomes are processed.
    """
    # Create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Define paths for temporary files inside the temp directory
        temp_bins = os.path.join(temp_dir, 'temp_bins.bed')
        temp_bins_filtered = os.path.join(temp_dir, 'temp_bins_filtered.bed')
        temp_bins_sorted = os.path.join(temp_dir, 'temp_bins_sorted.bed')
        temp_bedgraph_sorted = os.path.join(temp_dir, 'temp_bedgraph_sorted.bedgraph')
        temp_bedgraph_filtered = os.path.join(temp_dir, 'temp_bedgraph_filtered.bedgraph')
        temp_binned = os.path.join(temp_dir, 'temp_binned.bedgraph')

        # Use bedtools makewindows to create fixed-size bins
        with open(temp_bins, 'w') as f_bins:
            subprocess.run(['bedtools', 'makewindows', '-g', chrom_sizes, '-w', str(bin_size)],
                           stdout=f_bins, check=True)

        # Filter bins to specified chromosomes if chroms is provided
        if chroms:
            df_bins = pd.read_csv(temp_bins, sep='\t', header=None, names=['chrom', 'start', 'end'])
            df_bins = df_bins[df_bins['chrom'].isin(chroms)]
            if df_bins.empty:
                print(f"No bins found for the specified chromosomes: {chroms}")
                sys.exit(1)
            df_bins.to_csv(temp_bins_filtered, sep='\t', header=False, index=False)
        else:
            temp_bins_filtered = temp_bins

        # Sort temp_bins_filtered.bed
        with open(temp_bins_sorted, 'w') as f_bins_sorted:
            subprocess.run(['bedtools', 'sort', '-faidx', chrom_sizes, '-i', temp_bins_filtered],
                           stdout=f_bins_sorted, check=True)

        # Sort input_bedgraph
        with open(temp_bedgraph_sorted, 'w') as f_bedgraph_sorted:
            subprocess.run(['bedtools', 'sort', '-faidx', chrom_sizes, '-i', input_bedgraph],
                           stdout=f_bedgraph_sorted, check=True)

        # Filter input_bedgraph_sorted to specified chromosomes if chroms is provided
        if chroms:
            df_bedgraph = pd.read_csv(temp_bedgraph_sorted, sep='\t', header=None,
                                      names=['chrom', 'start', 'end', 'value'])
            df_bedgraph = df_bedgraph[df_bedgraph['chrom'].isin(chroms)]
            df_bedgraph.to_csv(temp_bedgraph_filtered, sep='\t', header=False, index=False)
            bedgraph_input = temp_bedgraph_filtered
        else:
            bedgraph_input = temp_bedgraph_sorted

        # Map the BedGraph data to the bins
        with open(temp_binned, 'w') as f_binned:
            subprocess.run(['bedtools', 'map', '-a', temp_bins_sorted, '-b', bedgraph_input,
                            '-c', '4', '-o', 'mean'], stdout=f_binned, check=True)

        # Read the binned data, specifying na_values
        df = pd.read_csv(
            temp_binned,
            sep='\t',
            header=None,
            names=['chrom', 'start', 'end', 'value'],
            na_values='.'
        )

        # Replace NaN values with zero
        df['value'] = df['value'].fillna(0)

        # Convert 'value' column to float
        df['value'] = df['value'].astype(float)

        # Save the binned BedGraph file
        df.to_csv(output_bedgraph, sep='\t', header=False, index=False)

        # The temporary directory and all its contents are deleted here



def localPearson_and_enrichmentSignificance(df, column1='readNum_1', column2='readNum_2', bin_number_of_window=11, step=1,
                                output_dir='output', chroms=None, **method_params):
    """
    Method to calculate local Pearson correlation and enrichment.
    This function processes the DataFrame and saves the output directly.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the merged tracks with columns 'chr', 'start', 'end', column1, column2.
    - column1 (str): Column name for track1 data in the DataFrame (default is 'readNum_1').
    - column2 (str): Column name for track2 data in the DataFrame (default is 'readNum_2').
    - bin_number_of_window (int): Number of bins in the sliding window. Default is 11.
    - step (int): Step size for the sliding window. Default is 1.
    - output_dir (str): Directory to save the output files. Default is 'output'.
    - chroms (list): List of chromosomes to process. If None, all chromosomes are processed.
    - **method_params: Additional method-specific parameters (e.g., 'percentile').

    Outputs:
    - Saves 'track_locCor.bedgraph' and 'track_ES.bedgraph' in the output directory.
    """
    print("Using localPearson_and_enrichmentSignificance (locP_and_ES)")
    # Extract method-specific parameters
    percentile = method_params.get('percentile', 5)  # Default percentile is 5

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file paths
    output_pearson_file = os.path.join(output_dir, 'track_locCor.bedgraph')
    output_wald_pvalue_file = os.path.join(output_dir, 'track_ES.bedgraph')

    # Initialize the final DataFrame
    df_final = pd.DataFrame()

    # Process specified chromosomes or all chromosomes
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = df['chr'].unique()

    for chrom in chromosomes:
        df_chr = df[df['chr'] == chrom].reset_index(drop=True)
        n = df_chr.shape[0]
        if n == 0:
            print(f"No data found for chromosome {chrom}")
            continue
        print(f"Processing {chrom} with {n} bins")

        # Replace values below the specified percentile within the chromosome
        all_values_chr = pd.concat([df_chr[column1], df_chr[column2]])
        non_zero_values_chr = all_values_chr[all_values_chr > 0]
        if len(non_zero_values_chr) == 0:
            percentile_value_chr = 0
        else:
            percentile_value_chr = np.percentile(non_zero_values_chr, percentile)
        print(f"{chrom} Percentile value: {percentile_value_chr}")

        df_chr[column1] = df_chr[column1].apply(lambda x: max(x, percentile_value_chr))
        df_chr[column2] = df_chr[column2].apply(lambda x: max(x, percentile_value_chr))

        # Ensure columns are in float64 to prevent the FutureWarning
        df_chr[column1] = df_chr[column1].astype('float64')
        df_chr[column2] = df_chr[column2].astype('float64')

        # Initialize lists to store calculated values
        mean_local_window_1 = [0] * n
        mean_local_window_2 = [0] * n
        var_local_window_1 = [0] * n
        var_local_window_2 = [0] * n
        dispersion_local_window_1 = [0] * n
        dispersion_local_window_2 = [0] * n
        pearson = [0] * n  # Initialize Pearson correlation values as floats
        df_chr_add = df_chr.copy()
        half_window = int((bin_number_of_window - 1) * 0.5)
        print("Half window size:", half_window)

        for i in range(half_window, n - half_window, step):
            window_1 = df_chr[column1].iloc[i - half_window: i + half_window + 1]
            window_2 = df_chr[column2].iloc[i - half_window: i + half_window + 1]
            mean_local_window_1[i] = float(window_1.mean())
            mean_local_window_2[i] = float(window_2.mean())
            var_local_window_1[i] = float(window_1.var(ddof=0))
            var_local_window_2[i] = float(window_2.var(ddof=0))
            dispersion_local_window_1[i] = max(
                (var_local_window_1[i] - mean_local_window_1[i]) / (mean_local_window_1[i] ** 2), 0)
            dispersion_local_window_2[i] = max(
                (var_local_window_2[i] - mean_local_window_2[i]) / (mean_local_window_2[i] ** 2), 0)
            if len(set(window_1)) > 1 and len(set(window_2)) > 1:
                pearson[i], _ = pearsonr(window_1, window_2)

        # Assign the results back to the DataFrame
        df_chr_add['mean_local_window_1'] = mean_local_window_1
        df_chr_add['mean_local_window_2'] = mean_local_window_2
        df_chr_add['var_local_window_1'] = var_local_window_1
        df_chr_add['var_local_window_2'] = var_local_window_2
        df_chr_add['pearson'] = pearson
        df_chr_add['dispersion_local_window_1'] = dispersion_local_window_1
        df_chr_add['dispersion_local_window_2'] = dispersion_local_window_2

        # Perform Wald test for the current chromosome
        n_chr = df_chr_add.shape[0]
        indices = df_chr_add.index[half_window:n_chr - half_window]

        # Initialize columns with float values to avoid dtype issues
        df_chr_add['SE'] = 0.0
        df_chr_add['log2Enrichment'] = 0.0
        df_chr_add['Wald'] = 0.0
        df_chr_add['Wald_pValue'] = 1.0
        df_chr_add['log_Wald_pValue'] = 0.0

        # Perform calculations with float64 values
        df_chr_add.loc[indices, 'SE'] = np.sqrt(
            1 / df_chr_add.loc[indices, 'mean_local_window_1'] + 1 / df_chr_add.loc[indices, 'mean_local_window_2'] +
            df_chr_add.loc[indices, 'dispersion_local_window_1'] + df_chr_add.loc[indices, 'dispersion_local_window_2']
        ).astype(float)

        df_chr_add.loc[indices, 'log2Enrichment'] = np.log2(
            df_chr_add.loc[indices, 'mean_local_window_2'] / df_chr_add.loc[indices, 'mean_local_window_1']
        ).astype(float)

        df_chr_add.loc[indices, 'Wald'] = df_chr_add.loc[indices, 'log2Enrichment'] / df_chr_add.loc[indices, 'SE']
        df_chr_add.loc[indices, 'Wald_pValue'] = 2 * (1 - norm.cdf(np.abs(df_chr_add.loc[indices, 'Wald'])))

        # Avoid division by zero or taking log of zero
        df_chr_add.loc[indices, 'log_Wald_pValue'] = -np.log10(
            df_chr_add.loc[indices, 'Wald_pValue'].replace(0, np.nan)
        ).fillna(0.0)

        # Append the processed chromosome data to df_final
        df_final = pd.concat([df_final, df_chr_add], ignore_index=True)

    if df_final.empty:
        print("No data processed.")
        sys.exit(1)

    # Update df with the processed data
    df = df_final

    # Extract the desired columns for output
    tracks_pearson = df[['chr', 'start', 'end', 'pearson']]
    tracks_log_Wald_pValue = df[['chr', 'start', 'end', 'log_Wald_pValue']]

    # Save the outputs
    tracks_pearson.to_csv(output_pearson_file, index=False, header=None, sep='\t')
    tracks_log_Wald_pValue.to_csv(output_wald_pvalue_file, index=False, header=None, sep='\t')

    print(f"Output files have been saved to {output_dir}")


def localWeightedPearson_and_enrichmentSignificance(df, column1='readNum_1', column2='readNum_2', bin_number_of_window=11, step=1,
                                        output_dir='output', chroms=None, **method_params):
    """
    Method to calculate local weighted Pearson correlation and enrichment.
    This function processes the DataFrame and saves the output directly.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the merged tracks with columns 'chr', 'start', 'end', column1, column2.
    - column1 (str): Column name for track1 data in the DataFrame (default is 'readNum_1').
    - column2 (str): Column name for track2 data in the DataFrame (default is 'readNum_2').
    - bin_number_of_window (int): Number of bins in the sliding window. Default is 11.
    - step (int): Step size for the sliding window. Default is 1.
    - output_dir (str): Directory to save the output files. Default is 'output'.
    - chroms (list): List of chromosomes to process. If None, all chromosomes are processed.
    - **method_params: Additional method-specific parameters (e.g., 'percentile').

    Outputs:
    - Saves 'track_locCor.bedgraph' and 'track_ES.bedgraph' in the output directory.
    """
    print("Using localWeightedPearson_and_enrichmentSignificance (locWP_and_ES)")
    # Extract method-specific parameters
    percentile = method_params.get('percentile', 5)  # Default percentile is 5

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file paths
    output_weightedPearson_file = os.path.join(output_dir, 'track_locCor.bedgraph')
    output_wald_pvalue_file = os.path.join(output_dir, 'track_ES.bedgraph')

    # Initialize the final DataFrame
    df_final = pd.DataFrame()

    # Process specified chromosomes or all chromosomes
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = df['chr'].unique()

    for chrom in chromosomes:
        df_chr = df[df['chr'] == chrom].reset_index(drop=True)
        n = df_chr.shape[0]
        if n == 0:
            print(f"No data found for chromosome {chrom}")
            continue
        print(f"Processing {chrom} with {n} bins")

        # Replace values below the specified percentile within the chromosome
        all_values_chr = pd.concat([df_chr[column1], df_chr[column2]])
        non_zero_values_chr = all_values_chr[all_values_chr > 0]
        if len(non_zero_values_chr) == 0:
            percentile_value_chr = 0
        else:
            percentile_value_chr = np.percentile(non_zero_values_chr, percentile)
        print(f"{chrom} Percentile value: {percentile_value_chr}")

        df_chr[column1] = df_chr[column1].apply(lambda x: max(x, percentile_value_chr))
        df_chr[column2] = df_chr[column2].apply(lambda x: max(x, percentile_value_chr))

        # Ensure columns are in float64 to prevent the FutureWarning
        df_chr[column1] = df_chr[column1].astype('float64')
        df_chr[column2] = df_chr[column2].astype('float64')

        mean_local_window_1 = [0] * n
        mean_local_window_2 = [0] * n
        var_local_window_1 = [0] * n
        var_local_window_2 = [0] * n
        dispersion_local_window_1 = [0] * n
        dispersion_local_window_2 = [0] * n
        weighted_pearson = [0] * n
        df_chr_add = df_chr.copy()
        half_window = int((bin_number_of_window - 1) * 0.5)
        print("Half window size:", half_window)
        weights_local = np.array([2 ** (-abs(i)) for i in range(-half_window, half_window + 1)])
        print("Weights_local array:", weights_local)

        for i in range(half_window, n - half_window, step):
            window_1 = df_chr[column1].iloc[i - half_window: i + half_window + 1]
            window_2 = df_chr[column2].iloc[i - half_window: i + half_window + 1]
            weighted_1 = window_1 * weights_local
            weighted_2 = window_2 * weights_local
            mean_local_window_1[i] = window_1.mean()
            mean_local_window_2[i] = window_2.mean()
            var_local_window_1[i] = window_1.var(ddof=0)
            var_local_window_2[i] = window_2.var(ddof=0)
            dispersion_local_window_1[i] = max(
                (var_local_window_1[i] - mean_local_window_1[i]) / (mean_local_window_1[i] ** 2), 0)
            dispersion_local_window_2[i] = max(
                (var_local_window_2[i] - mean_local_window_2[i]) / (mean_local_window_2[i] ** 2), 0)
            if len(set(window_1)) > 1 and len(set(window_2)) > 1:
                weighted_pearson[i], _ = pearsonr(weighted_1, weighted_2)
        df_chr_add['mean_local_window_1'] = mean_local_window_1
        df_chr_add['mean_local_window_2'] = mean_local_window_2
        df_chr_add['var_local_window_1'] = var_local_window_1
        df_chr_add['var_local_window_2'] = var_local_window_2
        df_chr_add['weighted_pearson'] = weighted_pearson
        df_chr_add['dispersion_local_window_1'] = dispersion_local_window_1
        df_chr_add['dispersion_local_window_2'] = dispersion_local_window_2

        # Perform Wald test for the current chromosome
        n_chr = df_chr_add.shape[0]
        indices = df_chr_add.index[half_window:n_chr - half_window]

        # Initialize columns
        df_chr_add['SE'] = 0.0
        df_chr_add['log2Enrichment'] = 0.0
        df_chr_add['Wald'] = 0.0
        df_chr_add['Wald_pValue'] = 1.0
        df_chr_add['log_Wald_pValue'] = 0.0

        # Perform calculations
        df_chr_add.loc[indices, 'SE'] = np.sqrt(
            1 / df_chr_add.loc[indices, 'mean_local_window_1'] + 1 / df_chr_add.loc[indices, 'mean_local_window_2'] +
            df_chr_add.loc[indices, 'dispersion_local_window_1'] + df_chr_add.loc[indices, 'dispersion_local_window_2']
        ).astype('float64')
        df_chr_add.loc[indices, 'log2Enrichment'] = np.log2(
            df_chr_add.loc[indices, 'mean_local_window_2'] / df_chr_add.loc[indices, 'mean_local_window_1']
        )
        df_chr_add.loc[indices, 'Wald'] = df_chr_add.loc[indices, 'log2Enrichment'] / df_chr_add.loc[indices, 'SE']
        df_chr_add.loc[indices, 'Wald_pValue'] = 2 * (1 - norm.cdf(np.abs(df_chr_add.loc[indices, 'Wald'])))

        # Avoid division by zero or taking log of zero
        df_chr_add.loc[indices, 'log_Wald_pValue'] = -np.log10(
            df_chr_add.loc[indices, 'Wald_pValue'].replace(0, np.nan)
        ).fillna(0)

        # Append the processed chromosome data to df_final
        df_final = pd.concat([df_final, df_chr_add], ignore_index=True)

    if df_final.empty:
        print("No data processed.")
        sys.exit(1)

    # Update df with the processed data
    df = df_final

    # Extract the desired columns for output
    tracks_weightedPearson = df[['chr', 'start', 'end', 'weighted_pearson']]
    tracks_log_Wald_pValue = df[['chr', 'start', 'end', 'log_Wald_pValue']]

    # Save the outputs
    tracks_weightedPearson.to_csv(output_weightedPearson_file, index=False, header=None, sep='\t')
    tracks_log_Wald_pValue.to_csv(output_wald_pvalue_file, index=False, header=None, sep='\t')

    print(f"Output files have been saved to {output_dir}")


def localSpearman_and_enrichmentSignificance(df, column1='readNum_1', column2='readNum_2', bin_number_of_window=11, step=1,
                                 output_dir='output', chroms=None, **method_params):
    """
    Method to calculate local Spearman correlation and enrichment.
    This function processes the DataFrame and saves the output directly.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the merged tracks with columns 'chr', 'start', 'end', column1, column2.
    - column1 (str): Column name for track1 data in the DataFrame (default is 'readNum_1').
    - column2 (str): Column name for track2 data in the DataFrame (default is 'readNum_2').
    - bin_number_of_window (int): Number of bins in the sliding window. Default is 11.
    - step (int): Step size for the sliding window. Default is 1.
    - output_dir (str): Directory to save the output files. Default is 'output'.
    - chroms (list): List of chromosomes to process. If None, all chromosomes are processed.
    - **method_params: Additional method-specific parameters (e.g., 'percentile').

    Outputs:
    - Saves 'track_locCor.bedgraph' and 'track_ES.bedgraph' in the output directory.
    """
    print("Using localSpearman_and_enrichmentSignificance (locS_and_ES)")
    # Extract method-specific parameters
    percentile = method_params.get('percentile', 5)  # Default percentile is 5

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file paths
    output_spearman_file = os.path.join(output_dir, 'track_locCor.bedgraph')
    output_wald_pvalue_file = os.path.join(output_dir, 'track_ES.bedgraph')

    # Initialize the final DataFrame
    df_final = pd.DataFrame()

    # Process specified chromosomes or all chromosomes
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = df['chr'].unique()

    for chrom in chromosomes:
        df_chr = df[df['chr'] == chrom].reset_index(drop=True)
        n = df_chr.shape[0]
        if n == 0:
            print(f"No data found for chromosome {chrom}")
            continue
        print(f"Processing {chrom} with {n} bins")

        # Replace values below the specified percentile within the chromosome
        all_values_chr = pd.concat([df_chr[column1], df_chr[column2]])
        non_zero_values_chr = all_values_chr[all_values_chr > 0]
        if len(non_zero_values_chr) == 0:
            percentile_value_chr = 0
        else:
            percentile_value_chr = np.percentile(non_zero_values_chr, percentile)
        print(f"{chrom} Percentile value: {percentile_value_chr}")

        df_chr[column1] = df_chr[column1].apply(lambda x: max(x, percentile_value_chr))
        df_chr[column2] = df_chr[column2].apply(lambda x: max(x, percentile_value_chr))

        # Ensure columns are in float64 to prevent the FutureWarning
        df_chr[column1] = df_chr[column1].astype('float64')
        df_chr[column2] = df_chr[column2].astype('float64')

        mean_local_window_1 = [0] * n
        mean_local_window_2 = [0] * n
        var_local_window_1 = [0] * n
        var_local_window_2 = [0] * n
        dispersion_local_window_1 = [0] * n
        dispersion_local_window_2 = [0] * n
        spearman = [0] * n
        df_chr_add = df_chr.copy()
        half_window = int((bin_number_of_window - 1) * 0.5)
        print("Half window size:", half_window)

        for i in range(half_window, n - half_window, step):
            window_1 = df_chr[column1].iloc[i - half_window: i + half_window + 1]
            window_2 = df_chr[column2].iloc[i - half_window: i + half_window + 1]
            mean_local_window_1[i] = window_1.mean()
            mean_local_window_2[i] = window_2.mean()
            var_local_window_1[i] = window_1.var(ddof=0)
            var_local_window_2[i] = window_2.var(ddof=0)
            dispersion_local_window_1[i] = max(
                (var_local_window_1[i] - mean_local_window_1[i]) / (mean_local_window_1[i] ** 2), 0)
            dispersion_local_window_2[i] = max(
                (var_local_window_2[i] - mean_local_window_2[i]) / (mean_local_window_2[i] ** 2), 0)
            if len(set(window_1)) > 1 and len(set(window_2)) > 1:
                spearman[i], _ = spearmanr(window_1, window_2)

        df_chr_add['mean_local_window_1'] = mean_local_window_1
        df_chr_add['mean_local_window_2'] = mean_local_window_2
        df_chr_add['var_local_window_1'] = var_local_window_1
        df_chr_add['var_local_window_2'] = var_local_window_2
        df_chr_add['spearman'] = spearman
        df_chr_add['dispersion_local_window_1'] = dispersion_local_window_1
        df_chr_add['dispersion_local_window_2'] = dispersion_local_window_2

        # Perform Wald test for the current chromosome
        n_chr = df_chr_add.shape[0]
        indices = df_chr_add.index[half_window:n_chr - half_window]

        # Initialize columns
        df_chr_add['SE'] = 0.0
        df_chr_add['log2Enrichment'] = 0.0
        df_chr_add['Wald'] = 0.0
        df_chr_add['Wald_pValue'] = 1.0
        df_chr_add['log_Wald_pValue'] = 0.0

        # Perform calculations
        df_chr_add.loc[indices, 'SE'] = np.sqrt(
            1 / df_chr_add.loc[indices, 'mean_local_window_1'] + 1 / df_chr_add.loc[indices, 'mean_local_window_2'] +
            df_chr_add.loc[indices, 'dispersion_local_window_1'] + df_chr_add.loc[indices, 'dispersion_local_window_2']
        )
        df_chr_add.loc[indices, 'log2Enrichment'] = np.log2(
            df_chr_add.loc[indices, 'mean_local_window_2'] / df_chr_add.loc[indices, 'mean_local_window_1']
        )
        df_chr_add.loc[indices, 'Wald'] = df_chr_add.loc[indices, 'log2Enrichment'] / df_chr_add.loc[indices, 'SE']
        df_chr_add.loc[indices, 'Wald_pValue'] = 2 * (1 - norm.cdf(np.abs(df_chr_add.loc[indices, 'Wald'])))

        # Avoid division by zero or taking log of zero
        df_chr_add.loc[indices, 'log_Wald_pValue'] = -np.log10(
            df_chr_add.loc[indices, 'Wald_pValue'].replace(0, np.nan)
        ).fillna(0)

        # Append the processed chromosome data to df_final
        df_final = pd.concat([df_final, df_chr_add], ignore_index=True)

    if df_final.empty:
        print("No data processed.")
        sys.exit(1)

    # Update df with the processed data
    df = df_final

    # Extract the desired columns for output
    tracks_spearman = df[['chr', 'start', 'end', 'spearman']]
    tracks_log_Wald_pValue = df[['chr', 'start', 'end', 'log_Wald_pValue']]

    # Save the outputs
    tracks_spearman.to_csv(output_spearman_file, index=False, header=None, sep='\t')
    tracks_log_Wald_pValue.to_csv(output_wald_pvalue_file, index=False, header=None, sep='\t')

    print(f"Output files have been saved to {output_dir}")


def localWeightedSpearman_and_enrichmentSignificance(df, column1='readNum_1', column2='readNum_2', bin_number_of_window=11, step=1,
                                         output_dir='output', chroms=None, **method_params):
    """
    Method to calculate local weighted Spearman correlation and enrichment.
    This function processes the DataFrame and saves the output directly.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the merged tracks with columns 'chr', 'start', 'end', column1, column2.
    - column1 (str): Column name for track1 data in the DataFrame (default is 'readNum_1').
    - column2 (str): Column name for track2 data in the DataFrame (default is 'readNum_2').
    - bin_number_of_window (int): Number of bins in the sliding window. Default is 11.
    - step (int): Step size for the sliding window. Default is 1.
    - output_dir (str): Directory to save the output files. Default is 'output'.
    - chroms (list): List of chromosomes to process. If None, all chromosomes are processed.
    - **method_params: Additional method-specific parameters (e.g., 'percentile').

    Outputs:
    - Saves 'track_locCor.bedgraph' and 'track_ES.bedgraph' in the output directory.
    """
    print("Using localWeightedSpearman_and_enrichmentSignificance (locWS_and_ES)")
    # Extract method-specific parameters
    percentile = method_params.get('percentile', 5)  # Default percentile is 5

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file paths
    output_weightedSpearman_file = os.path.join(output_dir, 'track_locCor.bedgraph')
    output_wald_pvalue_file = os.path.join(output_dir, 'track_ES.bedgraph')

    # Initialize the final DataFrame
    df_final = pd.DataFrame()

    # Process specified chromosomes or all chromosomes
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = df['chr'].unique()

    for chrom in chromosomes:
        df_chr = df[df['chr'] == chrom].reset_index(drop=True)
        n = df_chr.shape[0]
        if n == 0:
            print(f"No data found for chromosome {chrom}")
            continue
        print(f"Processing {chrom} with {n} bins")

        # Replace values below the specified percentile within the chromosome
        all_values_chr = pd.concat([df_chr[column1], df_chr[column2]])
        non_zero_values_chr = all_values_chr[all_values_chr > 0]
        if len(non_zero_values_chr) == 0:
            percentile_value_chr = 0
        else:
            percentile_value_chr = np.percentile(non_zero_values_chr, percentile)
        print(f"{chrom} Percentile value: {percentile_value_chr}")

        df_chr[column1] = df_chr[column1].apply(lambda x: max(x, percentile_value_chr))
        df_chr[column2] = df_chr[column2].apply(lambda x: max(x, percentile_value_chr))

        # Ensure columns are in float64 to prevent the FutureWarning
        df_chr[column1] = df_chr[column1].astype('float64')
        df_chr[column2] = df_chr[column2].astype('float64')

        mean_local_window_1 = [0] * n
        mean_local_window_2 = [0] * n
        var_local_window_1 = [0] * n
        var_local_window_2 = [0] * n
        dispersion_local_window_1 = [0] * n
        dispersion_local_window_2 = [0] * n
        weighted_spearman = [0] * n
        df_chr_add = df_chr.copy()
        half_window = int((bin_number_of_window - 1) * 0.5)
        print("Half window size:", half_window)
        weights_local = np.array([2 ** (-abs(i)) for i in range(-half_window, half_window + 1)])
        print("Weights_local array:", weights_local)

        for i in range(half_window, n - half_window, step):
            window_1 = df_chr[column1].iloc[i - half_window: i + half_window + 1]
            window_2 = df_chr[column2].iloc[i - half_window: i + half_window + 1]
            weighted_1 = window_1 * weights_local
            weighted_2 = window_2 * weights_local
            mean_local_window_1[i] = window_1.mean()
            mean_local_window_2[i] = window_2.mean()
            var_local_window_1[i] = window_1.var(ddof=0)
            var_local_window_2[i] = window_2.var(ddof=0)
            dispersion_local_window_1[i] = max(
                (var_local_window_1[i] - mean_local_window_1[i]) / (mean_local_window_1[i] ** 2), 0)
            dispersion_local_window_2[i] = max(
                (var_local_window_2[i] - mean_local_window_2[i]) / (mean_local_window_2[i] ** 2), 0)

            if len(set(window_1)) > 1 and len(set(window_2)) > 1:
                weighted_spearman[i], _ = spearmanr(weighted_1, weighted_2)

        df_chr_add['mean_local_window_1'] = mean_local_window_1
        df_chr_add['mean_local_window_2'] = mean_local_window_2
        df_chr_add['var_local_window_1'] = var_local_window_1
        df_chr_add['var_local_window_2'] = var_local_window_2
        df_chr_add['weighted_spearman'] = weighted_spearman
        df_chr_add['dispersion_local_window_1'] = dispersion_local_window_1
        df_chr_add['dispersion_local_window_2'] = dispersion_local_window_2

        # Perform Wald test for the current chromosome
        n_chr = df_chr_add.shape[0]
        indices = df_chr_add.index[half_window:n_chr - half_window]

        # Initialize columns
        df_chr_add['SE'] = 0.0
        df_chr_add['log2Enrichment'] = 0.0
        df_chr_add['Wald'] = 0.0
        df_chr_add['Wald_pValue'] = 1.0
        df_chr_add['log_Wald_pValue'] = 0.0

        # Perform calculations
        df_chr_add.loc[indices, 'SE'] = np.sqrt(
            1 / df_chr_add.loc[indices, 'mean_local_window_1'] + 1 / df_chr_add.loc[indices, 'mean_local_window_2'] +
            df_chr_add.loc[indices, 'dispersion_local_window_1'] + df_chr_add.loc[indices, 'dispersion_local_window_2']
        )
        df_chr_add.loc[indices, 'log2Enrichment'] = np.log2(
            df_chr_add.loc[indices, 'mean_local_window_2'] / df_chr_add.loc[indices, 'mean_local_window_1']
        )
        df_chr_add.loc[indices, 'Wald'] = df_chr_add.loc[indices, 'log2Enrichment'] / df_chr_add.loc[indices, 'SE']
        df_chr_add.loc[indices, 'Wald_pValue'] = 2 * (1 - norm.cdf(np.abs(df_chr_add.loc[indices, 'Wald'])))

        # Avoid division by zero or taking log of zero
        df_chr_add.loc[indices, 'log_Wald_pValue'] = -np.log10(
            df_chr_add.loc[indices, 'Wald_pValue'].replace(0, np.nan)
        ).fillna(0)

        # Append the processed chromosome data to df_final
        df_final = pd.concat([df_final, df_chr_add], ignore_index=True)

    if df_final.empty:
        print("No data processed.")
        sys.exit(1)

    # Update df with the processed data
    df = df_final

    # Extract the desired columns for output
    tracks_weightedSpearman = df[['chr', 'start', 'end', 'weighted_spearman']]
    tracks_log_Wald_pValue = df[['chr', 'start', 'end', 'log_Wald_pValue']]

    # Save the outputs
    tracks_weightedSpearman.to_csv(output_weightedSpearman_file, index=False, header=None, sep='\t')
    tracks_log_Wald_pValue.to_csv(output_wald_pvalue_file, index=False, header=None, sep='\t')

    print(f"Output files have been saved to {output_dir}")


def localMI_and_enrichmentSignificance(df, column1='readNum_1', column2='readNum_2', bin_number_of_window=11, step=1,
                           output_dir='output', chroms=None, **method_params):
    """
    Method to calculate local Mutual Information (MI) and enrichment.
    This function processes the DataFrame and saves the output directly.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the merged tracks with columns 'chr', 'start', 'end', column1, column2.
    - column1 (str): Column name for track1 data in the DataFrame (default is 'readNum_1').
    - column2 (str): Column name for track2 data in the DataFrame (default is 'readNum_2').
    - bin_number_of_window (int): Number of bins in the sliding window. Default is 11.
    - step (int): Step size for the sliding window. Default is 1.
    - output_dir (str): Directory to save the output files. Default is 'output'.
    - chroms (list): List of chromosomes to process. If None, all chromosomes are processed.
    - **method_params: Additional method-specific parameters (e.g., 'percentile').

    Outputs:
    - Saves 'track_locCor.bedgraph' and 'track_ES.bedgraph' in the output directory.
    """
    print("Using localMutualInformation_and_enrichmentSignificance (locMI_and_ES)")

    # Kraskov MI Function (defined inside the localMI_and_enrichmentSignificance function)
    def kraskov_mi(x, y, k=10):
        """
        Compute mutual information using the Kraskov method for nearest neighbors.

        Parameters:
        - x, y: Arrays of data points.
        - k: The number of nearest neighbors.
        """
        assert len(x) == len(y)
        N = len(x)
        x = np.array(x).reshape(N, 1)
        y = np.array(y).reshape(N, 1)
        data = np.hstack((x, y))

        # Compute distances to the k-th nearest neighbor in joint space
        tree = NearestNeighbors(metric='chebyshev')
        tree.fit(data)
        distances, _ = tree.kneighbors(n_neighbors=k + 1)
        eps = distances[:, k]  # Distance to k-th nearest neighbor

        # Initialize counts
        nx = np.zeros(N)
        ny = np.zeros(N)

        # Compute marginal counts
        tree_x = NearestNeighbors(metric='chebyshev')
        tree_x.fit(x)
        tree_y = NearestNeighbors(metric='chebyshev')
        tree_y.fit(y)

        for i in range(N):
            # Adjust epsilon to exclude points at distance exactly equal to eps[i]
            if eps[i] > 0:
                eps_i = eps[i] - 1e-10
            else:
                eps_i = 0.0

            # Count the number of neighbors within eps_i in the marginal spaces
            nx[i] = len(tree_x.radius_neighbors([x[i]], radius=eps_i, return_distance=False)[0]) - 1
            ny[i] = len(tree_y.radius_neighbors([y[i]], radius=eps_i, return_distance=False)[0]) - 1

        # Compute mutual information
        mi = digamma(k) - (np.mean(digamma(nx + 1) + digamma(ny + 1))) + digamma(N)
        return mi

    # Extract method-specific parameters
    percentile = method_params.get('percentile', 5)  # Default percentile is 5

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define the output file paths
    output_mutualInformation_file = os.path.join(output_dir, 'track_locCor.bedgraph')
    output_wald_pvalue_file = os.path.join(output_dir, 'track_ES.bedgraph')

    # Initialize the final DataFrame
    df_final = pd.DataFrame()

    # Process specified chromosomes or all chromosomes
    if chroms:
        chromosomes = chroms
    else:
        chromosomes = df['chr'].unique()

    for chrom in chromosomes:
        df_chr = df[df['chr'] == chrom].reset_index(drop=True)
        n = df_chr.shape[0]
        if n == 0:
            print(f"No data found for chromosome {chrom}")
            continue
        print(f"Processing {chrom} with {n} bins")

        # Replace values below the specified percentile within the chromosome
        all_values_chr = pd.concat([df_chr[column1], df_chr[column2]])
        non_zero_values_chr = all_values_chr[all_values_chr > 0]
        if len(non_zero_values_chr) == 0:
            percentile_value_chr = 0
        else:
            percentile_value_chr = np.percentile(non_zero_values_chr, percentile)
        print(f"{chrom} Percentile value: {percentile_value_chr}")

        df_chr[column1] = df_chr[column1].apply(lambda x: max(x, percentile_value_chr))
        df_chr[column2] = df_chr[column2].apply(lambda x: max(x, percentile_value_chr))

        # Ensure columns are in float64 to prevent the FutureWarning
        df_chr[column1] = df_chr[column1].astype('float64')
        df_chr[column2] = df_chr[column2].astype('float64')

        # Initialize variables
        mean_local_window_1 = [0] * n
        mean_local_window_2 = [0] * n
        var_local_window_1 = [0] * n
        var_local_window_2 = [0] * n
        dispersion_local_window_1 = [0] * n
        dispersion_local_window_2 = [0] * n
        mutual_information = [0] * n
        df_chr_add = df_chr.copy()
        half_window = int((bin_number_of_window - 1) * 0.5)
        print("Half window size:", half_window)

        for i in range(half_window, n - half_window, step):
            window_1 = df_chr[column1].iloc[i - half_window: i + half_window + 1]
            window_2 = df_chr[column2].iloc[i - half_window: i + half_window + 1]

            # Calculate Mutual Information
            mutual_information[i] = max(kraskov_mi(window_1, window_2, k=bin_number_of_window - 2), 0)

            # Compute statistics
            mean_local_window_1[i] = window_1.mean()
            mean_local_window_2[i] = window_2.mean()
            var_local_window_1[i] = window_1.var(ddof=0)
            var_local_window_2[i] = window_2.var(ddof=0)
            dispersion_local_window_1[i] = max(
                (var_local_window_1[i] - mean_local_window_1[i]) / (mean_local_window_1[i] ** 2), 0)
            dispersion_local_window_2[i] = max(
                (var_local_window_2[i] - mean_local_window_2[i]) / (mean_local_window_2[i] ** 2), 0)

        # Add results to the DataFrame
        df_chr_add['mean_local_window_1'] = mean_local_window_1
        df_chr_add['mean_local_window_2'] = mean_local_window_2
        df_chr_add['var_local_window_1'] = var_local_window_1
        df_chr_add['var_local_window_2'] = var_local_window_2
        df_chr_add['mutual_information'] = mutual_information
        df_chr_add['dispersion_local_window_1'] = dispersion_local_window_1
        df_chr_add['dispersion_local_window_2'] = dispersion_local_window_2

        # Perform Wald test for the current chromosome
        n_chr = df_chr_add.shape[0]
        indices = df_chr_add.index[half_window:n_chr - half_window]

        # Initialize columns
        df_chr_add['SE'] = 0.0
        df_chr_add['log2Enrichment'] = 0.0
        df_chr_add['Wald'] = 0.0
        df_chr_add['Wald_pValue'] = 1.0
        df_chr_add['log_Wald_pValue'] = 0.0

        # Perform calculations
        df_chr_add.loc[indices, 'SE'] = np.sqrt(
            1 / df_chr_add.loc[indices, 'mean_local_window_1'] + 1 / df_chr_add.loc[indices, 'mean_local_window_2'] +
            df_chr_add.loc[indices, 'dispersion_local_window_1'] + df_chr_add.loc[indices, 'dispersion_local_window_2']
        )
        df_chr_add.loc[indices, 'log2Enrichment'] = np.log2(
            df_chr_add.loc[indices, 'mean_local_window_2'] / df_chr_add.loc[indices, 'mean_local_window_1']
        )
        df_chr_add.loc[indices, 'Wald'] = df_chr_add.loc[indices, 'log2Enrichment'] / df_chr_add.loc[indices, 'SE']
        df_chr_add.loc[indices, 'Wald_pValue'] = 2 * (1 - norm.cdf(np.abs(df_chr_add.loc[indices, 'Wald'])))

        # Avoid division by zero or taking log of zero
        df_chr_add.loc[indices, 'log_Wald_pValue'] = -np.log10(
            df_chr_add.loc[indices, 'Wald_pValue'].replace(0, np.nan)
        ).fillna(0)

        df_final = pd.concat([df_final, df_chr_add])

    if df_final.empty:
        print("No data processed.")
        sys.exit(1)

    # Update df with the processed data
    df = df_final

    # Extract the desired columns for output
    tracks_mutualInformation = df[['chr', 'start', 'end', 'mutual_information']]
    tracks_log_Wald_pValue = df[['chr', 'start', 'end', 'log_Wald_pValue']]

    # Save the outputs
    tracks_mutualInformation.to_csv(output_mutualInformation_file, index=False, header=None, sep='\t')
    tracks_log_Wald_pValue.to_csv(output_wald_pvalue_file, index=False, header=None, sep='\t')

    print(f"Processing complete. Outputs saved to {output_dir}.")


def visualize_tracks(input_files, output_file, method='pyGenomeTracks', region=None, colors=None):
    if method == 'pyGenomeTracks':
        visualize_with_pygenometracks(
            input_files=input_files,
            output_file=output_file,
            region=region,
            colors=colors
        )
    elif method == 'plotly':
        visualize_with_plotly(
            input_files=input_files,
            output_file=output_file,
            region=region,
            colors=colors
        )
    else:
        raise ValueError(f"Unsupported visualization method: {method}")

def visualize_with_pygenometracks(input_files, output_file, region=None, colors=None):
    with tempfile.TemporaryDirectory() as temp_dir:
        config_file = os.path.join(temp_dir, 'tracks.ini')
        with open(config_file, 'w') as config:
            for idx, bedgraph_file in enumerate(input_files):
                track_id = f"track{idx + 1}"
                config.write(f'[{track_id}]\n')
                config.write(f'file = {bedgraph_file}\n')
                config.write(f'title = {os.path.basename(bedgraph_file)}\n')
                config.write('height = 4\n')
                config.write('type = line\n')

                # Use the provided or default color for each track
                try:
                    track_color = colors[idx]
                except (TypeError, IndexError):
                    print(f"Error: Insufficient colors provided for the number of input files.")
                    track_color = get_plotly_default_colors(1)[0]  # Fallback to default color

                config.write(f'color = {track_color}\n')
                config.write('\n')

        cmd = [
            'pyGenomeTracks',
            '--tracks', config_file,
            '--outFileName', output_file
        ]

        if region:
            chrom, start, end = region
            region_str = f'{chrom}:{start}-{end}'
            cmd.extend(['--region', region_str])

        try:
            subprocess.run(cmd, check=True)
            print(f"Visualization saved to {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running pyGenomeTracks: {e}")
            raise

def visualize_with_plotly(input_files, output_file, region=None, colors=None):
 # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
     
    num_tracks = len(input_files)
    fig = make_subplots(rows=num_tracks, cols=1, shared_xaxes=True, vertical_spacing=0.02)

    dataframes = []  # List to hold dataframes if needed for y-axis scaling

    for idx, bedgraph_file in enumerate(input_files):
        df = pd.read_csv(bedgraph_file, sep='\t', header=None)
        df.columns = ['chrom', 'start', 'end', 'value']

        # Filter by region if specified
        if region:
            chrom, start, end = region
            df = df[(df['chrom'] == chrom) & (df['end'] > start) & (df['start'] < end)]

        if df.empty:
            print(f"No data to plot for track {bedgraph_file}")
            continue

        df['position'] = (df['start'] + df['end']) / 2
        dataframes.append(df)

        try:
            track_color = colors[idx]
        except (TypeError, IndexError):
            print(f"Error: Insufficient colors provided for the number of input files.")
            track_color = get_plotly_default_colors(1)[0]  # Fallback to default color

        fig.add_trace(
            go.Scattergl(
                x=df['position'],
                y=df['value'],
                mode='lines',
                name=os.path.basename(bedgraph_file),
                line=dict(width=1, color=track_color)  # Apply color here
            ),
            row=idx + 1,
            col=1
        )

        # Add y-axis title for each subplot
        fig.update_yaxes(title_text=os.path.basename(bedgraph_file), row=idx + 1, col=1)

    fig.update_layout(
        title='Genomic Track Visualization',
        xaxis_title='Genomic Position',
        showlegend=False,
        height=300 * num_tracks,
        hovermode='x unified',
        plot_bgcolor='white',
    )

    # Update x-axis line color and thickness
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    # Update y-axis line color and thickness for all subplots
    for idx in range(num_tracks):
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black', row=idx + 1, col=1)

    fig.write_html(output_file)
    print(f"Visualization saved to {output_file}")

# Function to get the default Plotly color sequence for N tracks
def get_plotly_default_colors(num_colors):
    # Use Plotly's default color sequence (which has at least 10 colors by default)
    plotly_colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
        '#bcbd22', '#17becf'
    ]
    # Extend the color list if needed
    if num_colors > len(plotly_colors):
        # Repeat the color sequence to match the number of input files
        times = num_colors // len(plotly_colors) + 1
        plotly_colors = (plotly_colors * times)[:num_colors]
    return plotly_colors[:num_colors]


def find_significantly_different_regions(
        track_enrichment_file,
        track_correlation_file,
        output_dir,
        min_regionSize,
        enrichment_high_percentile,
        enrichment_low_percentile,
        corr_high_percentile,
        corr_low_percentile,
        chroms=None,
        chrom_sizes=None
):
    df_enrichment = pd.read_csv(
        track_enrichment_file, sep='\t', header=None, names=['chr', 'start', 'end', 'ES']
    )
    df_corr = pd.read_csv(
        track_correlation_file, sep='\t', header=None, names=['chr', 'start', 'end', 'locCorrelation']
    )

    # Merge the dataframes on genomic coordinates
    df = pd.merge(df_enrichment, df_corr, on=['chr', 'start', 'end'])

    # Remove bins with NaN values
    df.dropna(subset=['ES', 'locCorrelation'], inplace=True)

    # Determine chromosomes from chrom_sizes file if chroms='all' or chroms=None
    if chroms and chroms != ['all']:
        df = df[df['chr'].isin(chroms)]
        if df.empty:
            print(f"No data found for the specified chromosomes: {chroms}")
            sys.exit(1)
    else:
        # Extract chromosomes from the chrom_sizes file
        chromosomes = get_chromosomes_from_chrom_sizes(chrom_sizes)
        if not chromosomes:
            print("No chromosomes found in the chromosome sizes file.")
            sys.exit(1)
        df = df[df['chr'].isin(chromosomes)]
        if df.empty:
            print("No data found for chromosomes from chrom_sizes.")
            sys.exit(1)

    # Compute percentile ranks
    df['enrichment_percentile'] = df['ES'].rank(pct=True) * 100
    df['locCorrelation_percentile'] = df['locCorrelation'].rank(pct=True) * 100

    # Call visualize_scatter_ES_and_locCorr here
    visualize_scatter_ES_and_locCorr(df, output_dir, enrichment_high_percentile, enrichment_low_percentile,
                                     corr_high_percentile, corr_low_percentile)

    # Define thresholds
    # For enrichment
    df['high_enrichment'] = df['enrichment_percentile'] >= enrichment_high_percentile
    df['low_enrichment'] = df['enrichment_percentile'] <= enrichment_low_percentile

    # For correlation
    df['high_corr'] = df['locCorrelation_percentile'] >= corr_high_percentile
    df['low_corr'] = df['locCorrelation_percentile'] <= corr_low_percentile

    # Identify categories
    df['high_enrichment_high_corr'] = df['high_enrichment'] & df['high_corr']
    df['high_enrichment_low_corr'] = df['high_enrichment'] & df['low_corr']
    df['low_enrichment_high_corr'] = df['low_enrichment'] & df['high_corr']
    df['low_enrichment_low_corr'] = df['low_enrichment'] & df['low_corr']

    # Process each category
    categories = {
        'high_enrichment_high_correlation': df[df['high_enrichment_high_corr']],
        'high_enrichment_low_correlation': df[df['high_enrichment_low_corr']],
        'low_enrichment_high_correlation': df[df['low_enrichment_high_corr']],
        'low_enrichment_low_correlation': df[df['low_enrichment_low_corr']],
    }

    for category, df_category in categories.items():
        if df_category.empty:
            print(f"No regions found for category {category}.")
            continue

        # Save the individual significant bins for this category
        bins_output_file = os.path.join(output_dir, f"{category}_bins.bed")
        df_category[['chr', 'start', 'end', 'ES', 'locCorrelation']].to_csv(
            bins_output_file, sep='\t', header=True, index=False
        )
        print(f"Significant bins for category {category} saved to {bins_output_file}.")

        # Merge adjacent bins into regions
        significant_regions = merge_bins_into_regions(df_category, min_regionSize)

        if significant_regions.empty:
            print(f"No regions meet the minimum size for category {category}.")
            continue

        filename = f"{category}_regions.bed"
        output_path = os.path.join(output_dir, filename)
        significant_regions.to_csv(output_path, sep='\t', header=False, index=False)
        print(f"Found {len(significant_regions)} regions for category {category}.")
        print(f"Regions saved to {output_path}.")


def merge_bins_into_regions(df_bins, min_regionSize):
    """
    Merges adjacent bins into regions based on chromosome and start-end positions.

    Parameters:
    - df_bins (pd.DataFrame): DataFrame containing the bins to merge.
    - min_regionSize (int): Minimum number of consecutive bins to define a region.

    Returns:
    - df_regions (pd.DataFrame): DataFrame containing merged regions.
    """
    # Ensure the bins are sorted
    df_bins = df_bins.sort_values(['chr', 'start']).reset_index(drop=True)

    significant_regions = []
    current_region = None

    for idx, row in df_bins.iterrows():
        chrom = row['chr']
        start = row['start']
        end = row['end']

        if current_region is None:
            # Start a new region
            current_region = {'chr': chrom, 'start': start, 'end': end, 'bins': 1}
        else:
            # Check if the current bin is adjacent to the previous bin
            if chrom == current_region['chr'] and start == current_region['end']:
                # Extend the current region
                current_region['end'] = end
                current_region['bins'] += 1
            else:
                # Save the current region if it meets the minimum size
                if current_region['bins'] >= min_regionSize:
                    significant_regions.append(current_region)
                # Start a new region
                current_region = {'chr': chrom, 'start': start, 'end': end, 'bins': 1}

    # Add the last region if it meets the minimum size
    if current_region and current_region['bins'] >= min_regionSize:
        significant_regions.append(current_region)

    # Handle the case when significant_regions is empty
    if significant_regions:
        df_regions = pd.DataFrame(significant_regions)
    else:
        # Create an empty DataFrame with the required columns
        df_regions = pd.DataFrame(columns=['chr', 'start', 'end', 'bins'])

    return df_regions[['chr', 'start', 'end']]


def visualize_scatter_ES_and_locCorr(df, output_dir, enrichment_high_percentile, enrichment_low_percentile, corr_high_percentile,
                     corr_low_percentile):
    # Compute the actual threshold values
    enrichment_high_value = np.percentile(df['ES'], enrichment_high_percentile)
    enrichment_low_value = np.percentile(df['ES'], enrichment_low_percentile)
    corr_high_value = np.percentile(df['locCorrelation'], corr_high_percentile)
    corr_low_value = np.percentile(df['locCorrelation'], corr_low_percentile)

    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 10))
    gs = GridSpec(4, 4, figure=fig)

    # Main scatter plot
    ax_main = fig.add_subplot(gs[1:4, 0:3])
    scatter = ax_main.scatter(
        df['ES'],
        df['locCorrelation'],
        c='blue',
        alpha=0.5,
        edgecolor='k'
    )
    ax_main.set_xlabel('Enrich Significance')
    ax_main.set_ylabel('Local Correlation')
    ax_main.set_title('Enrich Significance vs Local Correlation')
    ax_main.grid(True)

    # Add threshold lines to scatter plot
    ax_main.axvline(
        x=enrichment_high_value,
        color='red',
        linestyle='--',
        label=f'Enrichment High Threshold ({enrichment_high_percentile}th Percentile)'
    )
    ax_main.axvline(
        x=enrichment_low_value,
        color='purple',
        linestyle='--',
        label=f'Enrichment Low Threshold ({enrichment_low_percentile}th Percentile)'
    )
    ax_main.axhline(
        y=corr_high_value,
        color='green',
        linestyle='--',
        label=f'Correlation High Threshold ({corr_high_percentile}th Percentile)'
    )
    ax_main.axhline(
        y=corr_low_value,
        color='orange',
        linestyle='--',
        label=f'Correlation Low Threshold ({corr_low_percentile}th Percentile)'
    )

    ax_main.legend(loc='upper right')

    # Cumulative distribution for Enrichment (ES) on top
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    sorted_enrichment = np.sort(df['ES'])
    cumulative_enrichment = np.arange(1, len(sorted_enrichment) + 1) / len(sorted_enrichment)
    ax_top.plot(sorted_enrichment, cumulative_enrichment, color='blue')
    ax_top.set_ylabel('Cumulative Distribution')
    ax_top.grid(True)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # Cumulative distribution for locCorrelation Correlation on the right
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    sorted_corr = np.sort(df['locCorrelation'])
    cumulative_corr = np.arange(1, len(sorted_corr) + 1) / len(sorted_corr)
    ax_right.plot(cumulative_corr, sorted_corr, color='blue')
    ax_right.set_xlabel('Cumulative Distribution')
    ax_right.grid(True)
    plt.setp(ax_right.get_yticklabels(), visible=False)

    # Adjust layout
    plt.subplots_adjust(wspace=0.05, hspace=0.05)

    # Save the plot
    plot_path = os.path.join(output_dir, 'ES_vs_locCorr.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Enrichment Significance vs Local Correlation scatter plot saved to {plot_path}.")
