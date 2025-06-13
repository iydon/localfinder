# File: utils.py
import tempfile, shutil, subprocess, sys, os, numpy as np, pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr, norm, spearmanr


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



def locCor_and_ES(df, column1='readNum_1', column2='readNum_2',
        bin_number_of_window=11, step=1, percentile=5, FC_thresh=1.5,
        bin_number_of_peak=11, corr_method='pearson',
        output_dir='output', chroms=None):

    """
    Calculate weighted local correlation and enrichment significance.
    Writes two bedgraph files: track_hmwC.bedgraph, track_ES.bedgraph.
    """

    print(f"calculate weighted local correlation and enrichment significance "
          f"(corr_method={corr_method})")
    print(f"parameters: percentile={percentile}, FC_thresh={FC_thresh}, "
          f"bin_number_of_window={bin_number_of_window}, "
          f"bin_number_of_peak={bin_number_of_peak}")
    EPS = 1e-9
    os.makedirs(output_dir, exist_ok=True)

    # ---------- 1  global normalisation ---------------------------------
    print("step1: global normalisation")
    cov1, cov2 = df[column1].sum(), df[column2].sum()
    if cov1 and cov2:
        if cov1 > cov2:
            df[column1] *= cov2 / cov1
            print(f"Scaled {column1} down to match {column2}")
        else:
            df[column2] *= cov1 / cov2
            print(f"Scaled {column2} down to match {column1}")

    # ---------- output paths -------------------------------------------
    out_ES   = os.path.join(output_dir, 'track_ES.bedgraph')
    out_hmwC = os.path.join(output_dir, 'track_hmwC.bedgraph')

    half_w = (bin_number_of_window - 1) // 2
    half_p = (bin_number_of_peak   - 1) // 2
    chromosomes = chroms if chroms else df['chr'].unique()

    df_final = pd.DataFrame()

    # ===================================================================
    # per-chromosome processing
    # ===================================================================
    print("step2: calculate weighted correlation and ES per chromosome")
    for chrom in chromosomes:
        df_raw = df[df['chr'] == chrom].reset_index(drop=True)           # raw (no floor)
        n = len(df_raw)
        if n == 0:
            print(f"No data for {chrom}, skipping")
            continue

        print(f"\nProcessing {chrom}: {n} bins")
        print(f"  half_window = {half_w}, half_peak = {half_p}")

        # 2·1  percentile floor (only for ES)
        all_vals = pd.concat([df_raw[column1], df_raw[column2]])
        nz = all_vals[all_vals > 0]
        pct = np.percentile(nz, percentile) if len(nz) else 0
        print(f"  {percentile}th percentile in non-zero bins = {pct:.4f}")

        df_floor = df_raw.copy()                                         ### floor
        mask = (df_floor[column1] <= pct) & (df_floor[column2] <= pct)
        df_floor.loc[mask, [column1, column2]] = pct 
        df_floor[column1] = df_floor[column1].astype('float64')
        df_floor[column2] = df_floor[column2].astype('float64')

        # 2·2  sliding-window statistics  (use *raw* data!)
        corr    = np.zeros(n)
        m_corr = np.zeros(n)
        hmw = np.zeros(n)

        m1_raw  = np.zeros(n); m2_raw  = np.zeros(n)    ### <<< CHANGED

        # floored parameters for logFC in peak region
        mES1_fl = np.zeros(n); mES2_fl = np.zeros(n)
        # floored parameters for SE in window region
        m1_fl = np.zeros(n);  m2_fl = np.zeros(n)                   ### <<< CHANGED
        v1_fl = np.zeros(n);  v2_fl = np.zeros(n)                   ### <<< CHANGED
        d1_fl = np.zeros(n);  d2_fl = np.zeros(n)                   ### <<< CHANGED

        print(f"Step2.2: calculate weighted local correlation ({corr_method})"
              f" and ES (window={bin_number_of_window}, step={step})")
        for i in range(half_w, n-half_w, step):
            # ---------- correlation window ----------  (RAW)  ###
            w1 = df_raw[column1].iloc[i-half_w : i+half_w+1]            ### <<< CHANGED (no-floor)
            w2 = df_raw[column2].iloc[i-half_w : i+half_w+1]            ### <<< CHANGED (no-floor)

            if len(set(w1)) > 1 and len(set(w2)) > 1:
                if corr_method == 'pearson':
                    r, _ = pearsonr(w1, w2)
                    corr[i] = max(r, 0)
                elif corr_method == 'spearman':
                    r, _ = spearmanr(w1, w2)
                    corr[i] = max(r, 0)
                else:
                    raise ValueError(f"Unknown corr_method: {corr_method}")

            m1_raw[i], m2_raw[i] = w1.mean(), w2.mean()
            hmw[i]  = (m1_raw[i]*m2_raw[i])/(m1_raw[i]+m2_raw[i]+EPS)

            # SE in window
            wf1 = df_floor[column1].iloc[i-half_w : i+half_w+1]
            wf2 = df_floor[column2].iloc[i-half_w : i+half_w+1]

            m1_fl[i], m2_fl[i] = wf1.mean(), wf2.mean()               ### <<< CHANGED
            v1_fl[i], v2_fl[i] = wf1.var(ddof=0), wf2.var(ddof=0)     ### <<< CHANGED
            d1_fl[i] = max((v1_fl[i] - m1_fl[i])/(m1_fl[i]**2 + EPS), 0)  ### <<< CHANGED
            d2_fl[i] = max((v2_fl[i] - m2_fl[i])/(m2_fl[i]**2 + EPS), 0)  ### <<< CHANGED
          

            # ---------- logFC in peak region ----------  (FLOORED)  ###
            p1 = df_floor[column1].iloc[i-half_p : i+half_p+1]          ### <<< CHANGED (no-floor)
            p2 = df_floor[column2].iloc[i-half_p : i+half_p+1]          ### <<< CHANGED (no-floor)
            mES1_fl[i], mES2_fl[i] = p1.mean(), p2.mean()

        # 2·3  calculate mean correlation
        for i in range(half_w, n-half_w, step):
            local_c = corr[i-half_w : i+half_w+1]
            m_corr[i] = local_c.mean()


        # ---------- logFC in peak region -----------------------------
        safe1 = np.maximum(mES1_fl, pct)
        safe2 = np.maximum(mES2_fl, pct)
        logFC = np.log(safe2 / (safe1 + EPS)) / np.log(FC_thresh)
        logFC[:half_w]  = 0.0
        logFC[-half_w:] = 0.0

        idx = np.arange(half_w, n-half_w, step)
        mu1 = np.maximum(m1_fl[idx], pct)                 
        mu2 = np.maximum(m2_fl[idx], pct)
        SE  = np.sqrt(1/mu1 + 1/mu2 + d1_fl[idx] + d2_fl[idx])        ### <<< CHANGED
        Wald = logFC[idx] / SE
        p    = 2 * (1 - norm.cdf(np.abs(Wald)))
        lP   = -np.log10(np.where(p == 0, np.nan, p))

        # ---------- assemble dataframe ----------
        dfc = df_raw.copy()
        dfc['corr']        = corr
        dfc['m_corr']      = m_corr
        dfc['hmw']         = hmw
        dfc['logFC']       = logFC
        dfc['SE'] = 0.0
        dfc['Wald']            = 0.0
        dfc['Wald_pValue']     = 0.0
        dfc['log_Wald_pValue'] = 0.0
        dfc.loc[idx, 'SE']            = SE
        dfc.loc[idx, 'Wald']          = Wald
        dfc.loc[idx, 'Wald_pValue']   = p
        dfc.loc[idx, 'log_Wald_pValue'] = np.nan_to_num(lP, nan=0.0)

        df_final = pd.concat([df_final, dfc], ignore_index=True)

    # ---------- 3  stack & write  --------------------------------------
    print("step3: write track_hmwC, track_ES")
    df_final['signed_log_Wald_pValue'] = (np.sign(df_final['logFC']) * df_final['log_Wald_pValue'])
    df_final['hmwC']   = df_final['m_corr'] * df_final['hmw']

    df_final[['chr', 'start', 'end', 'hmwC']].to_csv(out_hmwC,sep='\t', header=False, index=False)
    df_final[['chr', 'start', 'end', 'signed_log_Wald_pValue']].to_csv(out_ES,sep='\t', header=False, index=False)
    
    print(f"Saved outputs to {output_dir}")



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
    track_ES_file: str,
    track_hmwC_file: str,
    output_dir: str,
    p_thresh: float = 0.05,
    binNum_thresh: int = 2,
    chroms=None,
    chrom_sizes=None
):
    """
    1) Read track_ES.bedgraph and track_hmwC.bedgraph  
    2) Merge on chr/start/end  
    3) Keep bins with abs(ES) >= -log10(p_thresh) AND same sign runs of length>=min_region  
    4) Write those bins (chr, start, end, ES, hmwC) to '<output_dir>/signif_bins.tsv'  
    5) Merge adjacent bins into regions (same chr, end==next start), average hmwC across region  
       and sort descending by avg_hmwC, write to '<output_dir>/signif_regions.tsv'  
    6) Plot rank vs log(hmwC+1) for those regions, save '<output_dir>/hmwC_rank.png'  
    """
    os.makedirs(output_dir, exist_ok=True)

    # 1) read  
    df_es = pd.read_csv(track_ES_file, sep='\t', header=None,
                        names=['chr','start','end','ES'])
    print(df_es.shape)
    # print(df_es.head())
    df_c  = pd.read_csv(track_hmwC_file, sep='\t', header=None,
                        names=['chr','start','end','hmwC'])
    print(df_c.shape)
    # print(df_c.head())
    # 2) merge  
    df = pd.merge(df_es, df_c, on=['chr','start','end'])
    print(df.shape)
    # print(df.head())

    # 3) filter by ES  
    thresh = -np.log10(p_thresh)
    print(f'thresh: {thresh}')  
    # split positive and negative extremes
    df_pos = df[df['ES'] >=  thresh]
    df_neg = df[df['ES'] <= -thresh]
    # print(df_pos.shape)
    # print(df_pos.head(50))
    # print(df_neg.shape)
    # print(df_neg.head(50))

    # helper to merge runs on same chromosome
    def _merge_runs(df_bins):
        regs = []
        for chrom, grp in df_bins.groupby('chr'):
            grp = grp.sort_values('start')
            curr = None
            for _, row in grp.iterrows():
                if curr is None or row.start != curr['end']:
                    # flush old
                    if curr and curr['count'] >= binNum_thresh:
                        regs.append(curr)
                    curr = dict(chr=chrom,
                                start=row.start,
                                end=row.end,
                                wvals=[row.hmwC],
                                count=1)
                else:
                    curr['end']   = row.end
                    curr['wvals'].append(row.hmwC)
                    curr['count'] += 1
            # last
            if curr and curr['count'] >= binNum_thresh:
                regs.append(curr)
        # build DataFrame
        df_regs = pd.DataFrame(regs)
        # if no runs passed the threshold, just return empty frame with the right columns
        if df_regs.empty:
            return pd.DataFrame(columns=['chr','start','end','avg_hmwC'])
        df_regs['avg_hmwC'] = df_regs['wvals'].apply(np.mean)
        return df_regs[['chr','start','end','avg_hmwC']]

    # 4) write bins
    df_pos_sig = df_pos.sort_values(['chr','start'])
    pos_sig_n = df_pos_sig.shape[0]
    pos_bins_out = os.path.join(output_dir, 'pos_signif_bins.tsv')
    df_pos_sig[['chr','start','end','ES','hmwC']].to_csv(
        pos_bins_out, sep='\t', index=False
    )
    print(f"{pos_sig_n} pos ignificant bins written to {pos_bins_out}")
    df_neg_sig = df_neg.sort_values(['chr','start'])
    neg_sig_n = df_neg_sig.shape[0]
    neg_bins_out = os.path.join(output_dir, 'neg_signif_bins.tsv')
    df_neg_sig[['chr','start','end', 'ES','hmwC']].to_csv(
        neg_bins_out, sep='\t', index=False
    )
    print(f"{neg_sig_n} neg ignificant bins written to {neg_bins_out}")

    # 5) merge regions & write
    df_pos_regs = _merge_runs(df_pos).sort_values('avg_hmwC', ascending=False)
    pos_n = df_pos_regs.shape[0]
    print(pos_n)
    print(df_pos_regs.shape)
    print(df_pos_regs)

    pos_regs_out = os.path.join(output_dir, 'signif_pos_regions.tsv')
    df_pos_regs.to_csv(pos_regs_out, sep='\t', index=False)
    print(f"{pos_n} merged pos regions written to {pos_regs_out}")


    df_neg_regs = _merge_runs(df_neg).sort_values('avg_hmwC', ascending=False)
    neg_n = df_neg_regs.shape[0]
    print(neg_n)
    print(df_neg_regs.shape)
    print(df_neg_regs)

    neg_regs_out = os.path.join(output_dir, 'signif_neg_regions.tsv')
    df_neg_regs.to_csv(neg_regs_out, sep='\t', index=False)
    print(f"{neg_n} merged neg regions written to {neg_regs_out}")

    # 6) plot rank vs log(w+1)
    df_pos_regs['rank'] = np.arange(1, len(df_pos_regs)+1)
    plt.figure(figsize=(6,4))
    plt.plot(df_pos_regs['rank'], np.log(df_pos_regs['avg_hmwC']+1), marker='.', linestyle='none',markersize=1, color='k')
    plt.xlabel('Region rank')
    plt.ylabel('log(hmwC+1)')
    plt.title('Ranked region hmwC')
    plt.tight_layout()
    fig_out = os.path.join(output_dir, 'hmwC_rank_in_pos_regions.png')
    plt.savefig(fig_out, dpi=600)
    plt.close()
    print(f"Rank plot of pos regions saved to {fig_out}")

    df_neg_regs['rank'] = np.arange(1, len(df_neg_regs)+1)
    plt.figure(figsize=(6,4))
    plt.plot(df_neg_regs['rank'], np.log(df_neg_regs['avg_hmwC']+1), marker='.', linestyle='none',markersize=1, color='k')
    plt.xlabel('Region rank')
    plt.ylabel('log(hmwC+1)')
    plt.title('Ranked region hmwC')
    plt.tight_layout()
    fig_out = os.path.join(output_dir, 'hmwC_rank_in_neg_regions.png')
    plt.savefig(fig_out, dpi=600)
    plt.close()
    print(f"Rank plot of neg regions saved to {fig_out}")


