# localfinder

localfinder – calculate weighted local correlation (hmC) and enrichment significance (ES) between two genomic tracks, optionally discover significantly different regions, and visualise results. 

## Installation Requirements

Before installing and using `localfinder`, please ensure that the following external tools are installed on your `$PATH`:

- **bedtools**: Used for genomic interval operations.
  - Installation: [https://bedtools.readthedocs.io/en/latest/content/installation.html](https://bedtools.readthedocs.io/en/latest/content/installation.html)
  - conda install -c bioconda -c conda-forge bedtools 
  - mamba install -c bioconda -c conda-forge bedtools
- **bigWigToBedGraph (UCSC utility)**: Used for converting BigWig files to BedGraph format.
  - Download: [http://hgdownload.soe.ucsc.edu/admin/exe/](http://hgdownload.soe.ucsc.edu/admin/exe/)
  - conda install -c bioconda -c conda-forge ucsc-bigwigtobedgraph
  - mamba install -c bioconda -c conda-forge ucsc-bigwigtobedgraph
- **samtools**: Used for processing SAM/BAM files.
  - Installation: [http://www.htslib.org/download/](http://www.htslib.org/download/)
  - conda install -c bioconda -c conda-forge samtools
  - mamba install -c bioconda -c conda-forge samtools

These tools are required for processing genomic data and must be installed separately.

## Installation

Install `localfinder` using `pip`:

```bash
pip install localfinder
```

## Usage
There are 5 subcommands (bin, calc, findreg, viz, pipeline) in localfinder, and you can check it using:
```bash
localfinder -h
```
### bin
```bash
localfinder bin -h
```
```
usage: localfinder bin [-h] --input_files INPUT_FILES [INPUT_FILES ...] --output_dir OUTPUT_DIR [--bin_size BIN_SIZE] --chrom_sizes CHROM_SIZES [--chroms CHROMS [CHROMS ...]] [--threads N]

Bin genomic tracks into fixed-size bins and output BedGraph format.

optional arguments:
  -h, --help            show this help message and exit
  --input_files INPUT_FILES [INPUT_FILES ...]
                        Input files in BigWig/BedGraph/BAM/SAM format.
  --output_dir OUTPUT_DIR
                        Output directory for binned data.
  --bin_size BIN_SIZE   Size of each bin (default: 200).
  --chrom_sizes CHROM_SIZES
                        Path to the chromosome sizes file.
  --chroms CHROMS [CHROMS ...]
                        'all' or specific chromosomes (e.g. chr1 chr2).
  --threads N, -t N     Number of worker processes to run in parallel (default: 1).

Usage Example 1:
    localfinder bin --input_files track1.bw track2.bw --output_dir ./binned_tracks --bin_size 200 --chrom_sizes mm10.chrom.sizes --chroms chr1 chr2

Usage Example 2:
    localfinder bin --input_files track1.bigwig track2.bigwig --output_dir ./binned_tracks --bin_size 200 --chrom_sizes hg19.chrom.sizes --chroms all --threads 4
```

### calc
```bash
localfinder calc -h
```
```
usage: localfinder calc [-h] --track1 TRACK1 --track2 TRACK2 --output_dir OUTPUT_DIR [--method {locP_and_ES,locS_and_ES}] [--FDR] [--binNum_window BINNUM_WINDOW] [--step STEP]
                        [--percentile PERCENTILE] [--percentile_mode {all,nonzero}] [--binNum_peak BINNUM_PEAK] [--FC_thresh FC_THRESH] [--norm_method {scale,cpm}] --chrom_sizes CHROM_SIZES
                        [--chroms CHROMS [CHROMS ...]] [--threads N]

Calculate weighted local correlation and enrichment significance between two BedGraph tracks.

optional arguments:
  -h, --help            show this help message and exit
  --track1 TRACK1       First input BedGraph file.
  --track2 TRACK2       Second input BedGraph file.
  --output_dir OUTPUT_DIR
                        Output directory for results.
  --method {locP_and_ES,locS_and_ES}
                        Method for calculate weighted local correlation and enrichment significance. P and S denote Pearson and Spearman respectively (default: locP_and_ES).
  --FDR                 Use Benjamini–Hochberg FDR-corrected q-values instead of raw P-values.
  --binNum_window BINNUM_WINDOW
                        Number of bins in the sliding window (default: 11).
  --step STEP           Step size for the sliding window (default: 1)
  --percentile PERCENTILE
                        Percentile for floor correction of low-coverage bins (default: 90). High percentile such as 95 or 98 is recommended, when tracks mainly contains some high sharp peaks,
                        while small percentile like 90 is recommended when tracks mainly contain broad and relatively low peaks.
  --percentile_mode {all,nonzero}
                        Use all bins or only non-zero bins for percentile (default: all). When choosing 'all' mode, we will choose larger percentile value than that of 'nonzero' mode to get
                        the same floor.
  --binNum_peak BINNUM_PEAK
                        Number of bins of the peak for ES (default: 3). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended.
  --FC_thresh FC_THRESH
                        Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to
                        the FC_thresh in the vocalno plot. Wald, a statistical value following a normal distribution here, euqal to log1.5(track1/track2) / SE can be used to calculate the p
                        value, whose -log10 represents for ES here.
  --norm_method {scale,cpm}
                        Normalisation: scale (match totals) or cpm (counts per million).
  --chrom_sizes CHROM_SIZES
                        Path to the chromosome sizes file.
  --chroms CHROMS [CHROMS ...]
                        Chromosomes to process (e.g., chr1 chr2). Defaults to 'all'.
  --threads N, -t N     Worker processes per-chromosome (default: 1).

Usage Example 1:
    localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --output_dir ./results --method locP_and_ES --FDR --binNum_window 11 --step 1 --percentile 90 --percentile_mode all --binNum_peak 3 --FC_thresh 1.5 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2 --threads 4

Usage Example 2:
    localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --output_dir ./results --percentile 99 --binNum_peak 2 --chrom_sizes hg19.chrom.sizes
```

### findreg
```bash
localfinder findreg -h
```
```
usage: localfinder findreg [-h] --track_E TRACK_E --track_C TRACK_C --output_dir OUTPUT_DIR [--p_thresh P_THRESH] [--binNum_thresh BINNUM_THRESH] [--chroms CHROMS [CHROMS ...]] --chrom_sizes
                           CHROM_SIZES

Merge consecutive significant bins into regions. And find significantly different regions from ES & hmwC tracks.

optional arguments:
  -h, --help            show this help message and exit
  --track_E TRACK_E     track_ES.bedgraph
  --track_C TRACK_C     track_hmwC.bedgraph
  --output_dir OUTPUT_DIR
  --p_thresh P_THRESH   P-value threshold (default: 0.05)
  --binNum_thresh BINNUM_THRESH
                        Min consecutive significant bins per region (default: 2)
  --chroms CHROMS [CHROMS ...]
  --chrom_sizes CHROM_SIZES

Example 1:
  localfinder findreg --track_E track_ES.bedgraph --track_C track_hmwC.bedgraph --output_dir ./findreg_out --p_thresh 0.05 --binNum_thresh 2 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

Example 2:
  localfinder findreg --track_E track_ES.bedgraph --track_C track_hmwC.bedgraph --output_dir ./findreg_out --chrom_sizes hg19.chrom.sizes
```
### pipeline
```bash
localfinder pipeline -h
```
```
usage: localfinder findreg [-h] --track_E TRACK_E --track_C TRACK_C --output_dir OUTPUT_DIR [--p_thresh P_THRESH] [--binNum_thresh BINNUM_THRESH] [--chroms CHROMS [CHROMS ...]] --chrom_sizes
                           CHROM_SIZES

Merge consecutive significant bins into regions. And find significantly different regions from ES & hmwC tracks.

optional arguments:
  -h, --help            show this help message and exit
  --track_E TRACK_E     track_ES.bedgraph
  --track_C TRACK_C     track_hmwC.bedgraph
  --output_dir OUTPUT_DIR
  --p_thresh P_THRESH   P-value threshold (default: 0.05)
  --binNum_thresh BINNUM_THRESH
                        Min consecutive significant bins per region (default: 2)
  --chroms CHROMS [CHROMS ...]
  --chrom_sizes CHROM_SIZES

Example 1:
  localfinder findreg --track_E track_ES.bedgraph --track_C track_hmwC.bedgraph --output_dir ./findreg_out --p_thresh 0.05 --binNum_thresh 2 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

Example 2:
  localfinder findreg --track_E track_ES.bedgraph --track_C track_hmwC.bedgraph --output_dir ./findreg_out --chrom_sizes hg19.chrom.sizes
(tools_localfinder_v20) fayevalentine@Mac localfinder %localfinder pipeline -h
usage: localfinder pipeline [-h] --input_files INPUT_FILES [INPUT_FILES ...] --output_dir OUTPUT_DIR --chrom_sizes CHROM_SIZES [--bin_size BIN_SIZE] [--method {locP_and_ES,locS_and_ES}]
                            [--FDR] [--binNum_window BINNUM_WINDOW] [--binNum_peak BINNUM_PEAK] [--step STEP] [--percentile PERCENTILE] [--percentile_mode {all,nonzero}]
                            [--FC_thresh FC_THRESH] [--norm_method {scale,cpm}] [--p_thresh P_THRESH] [--binNum_thresh BINNUM_THRESH] [--chroms CHROMS [CHROMS ...]] [--threads THREADS]

Run localfinder sequentially.

optional arguments:
  -h, --help            show this help message and exit
  --input_files INPUT_FILES [INPUT_FILES ...]
                        Input files in BigWig/BedGraph format.
  --output_dir OUTPUT_DIR
                        Output directory for the pipeline results.
  --chrom_sizes CHROM_SIZES
                        Path to the chromosome sizes file.
  --bin_size BIN_SIZE   Size of each bin (default: 200).
  --method {locP_and_ES,locS_and_ES}
                        Method for calculate weighted local correlation and enrichment significance. P and S denote Pearson and Spearman respectively (default: locP_and_ES).
  --FDR                 Use Benjamini–Hochberg FDR-corrected q-values instead of raw P-values
  --binNum_window BINNUM_WINDOW
                        Number of bins in the sliding window (default: 11).
  --binNum_peak BINNUM_PEAK
                        Number of bins of the peak for ES (default: 3). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended.
  --step STEP           Step size for the sliding window (default: 1)
  --percentile PERCENTILE
                        Percentile for floor correction of low-coverage bins (default: 90). High percentile such as 95 or 98 is recommended, when tracks mainly contains some high sharp peaks,
                        while small percentile like 90 is recommended when tracks mainly contain broad and relatively low peaks.
  --percentile_mode {all,nonzero}
                        Use all bins or only non-zero bins for percentile (default: all). When choosing "all" mode, we will choose larger percentile value than that of "nonzero" mode to get
                        the same floor.
  --FC_thresh FC_THRESH
                        Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to
                        the FC_thresh in the volcano plot. Wald, a statistical value following a normal distribution here, equal to log1.5(track1/track2) / SE can be used to calculate the p
                        value, whose -log10 represents for ES here.
  --norm_method {scale,cpm}
                        Normalisation: scale (match totals) or cpm (counts per million).
  --p_thresh P_THRESH   P-value threshold for merging significant bins into regions (default: 0.05)
  --binNum_thresh BINNUM_THRESH
                        Min consecutive significant bins per region
  --chroms CHROMS [CHROMS ...]
                        Chromosomes to process (e.g., chr1 chr2). Defaults to "all".
  --threads THREADS, -t THREADS
                        Number of worker processes for bin & calc (default: 1)

Usage Example 1:
    localfinder pipeline --input_files track1.bedgraph track2.bedgraph --output_dir ./results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --FDR --binNum_window 11 --binNum_peak 3 --step 1 --percentile 90 --percentile_mode all --FC_thresh 1.5 --norm_method cpm --p_thresh 0.05 --binNum_thresh 2 --chroms chr1 chr2 --threads 4

Usage Example 2:
    localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./results --chrom_sizes hg19.chrom.sizes --binNum_peak 3 --percentile 95 --binNum_thresh 2
```
### viz
```bash
localfinder viz -h
```
```
    usage: localfinder viz [-h] --input_files INPUT_FILES [INPUT_FILES ...] --output_file OUTPUT_FILE --method {pyGenomeTracks,plotly} [--region CHROM START END] [--colors COLORS [COLORS ...]]
    
    Visualize genomic tracks.
    
    options:  
    -h, --help                                          `show this help message and exit`  
    --input_files INPUT_FILES [INPUT_FILES ...]         `Input BedGraph files to visualize`  
    --output_file OUTPUT_FILE                           `Output visualization file (e.g., PNG, HTML)`  
    --method {pyGenomeTracks,plotly}                    `Visualization method to use`  
    --region CHROM START END                            `Region to visualize in the format: CHROM START END (e.g., chr20 1000000 2000000)`  
    --colors COLORS [COLORS ...]                        `Colors for the tracks (optional)`
    
    Usage Example 1:
        localfinder viz --input_files track1.bedgraph track2.bedgraph --output_file output.html --method plotly --region chr1 1000000 2000000 --colors blue red
    
    Usage Example 2:
        localfinder viz --input_files track1.bedgraph track2.bedgraph --output_file output.png --method pyGenomeTracks --region chr1 1000000 2000000
```

## Run an example step by step
Create a conda env called localfinder and enter this conda environment
```bash
conda create -n localfinder
conda activate  localfinder
```

Install external tools and localfinder
```bash
conda install -c conda-forge -c bioconda samtools bedtools ucsc-bigwigtobedgraph
pip install localfinder
```

Download the souce code of [localfinder](https://github.com/astudentfromsustech/localfinder)  
```bash
git clone https://github.com/astudentfromsustech/localfinder.git
```

Run the examples under localfinder/tests/ (scripts have been preprared in tests folder)  
