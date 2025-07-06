# localfinder

A tool calculating weighted local correlation and enrichment significance of two tracks and finding significantly different genomic regions.

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
  - Warning: On Apple Silicon (osx-arm64) Bioconda doesn’t yet build the UCSC “bigWigToBedGraph” utility, so conda/mamba can’t find ucsc-bigwigtobedgraph. localfinder/tests/tools_bigWigToBedGraph_forARM/install_bigWigToBedGraph_forARM64.sh is a shell script to install it.
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
    usage: localfinder bin [-h] --input_files INPUT_FILES [INPUT_FILES ...] --output_dir OUTPUT_DIR --chrom_sizes CHROM_SIZES [--bin_size BIN_SIZE] [--chroms CHROMS [CHROMS ...]]
    
    Bin genomic tracks into fixed-size bins and output BedGraph format.
    
    options:  
      -h, --help                                        `Show this help message and exit`  
      --input_files INPUT_FILES [INPUT_FILES ...]       `Input files in BigWig/BedGraph/BAM/SAM format`    
      --output_dir OUTPUT_DIR                           `Output directory`
      --chrom_sizes CHROM_SIZES                         `Path to the chromosome sizes file`   
      --bin_size BIN_SIZE                               `Size of each bin (default: 200)`  
      --chroms CHROMS [CHROMS ...]                      `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`
    
    Usage Example 1:
        localfinder bin --input_files track1.bw track2.bw --output_dir ./bin_out --chrom_sizes hg19.chrom.sizes --bin_size 200 --chroms chr1 chr2
    
    Usage Example 2:
        localfinder bin --input_files track1.bigwig track2.bigwig --output_dir ./bin_out --chrom_sizes hg19.chrom.sizes --bin_size 200 --chroms all


### calc
```bash
localfinder calc -h
```

    usage: localfinder calc [-h] --track1 TRACK1 --track2 TRACK2 [--method {locP_and_ES,locS_and_ES}] [--FDR] [--binNum_window BINNUM_WINDOW] [--step STEP]
                        [--percentile PERCENTILE] [--binNum_peak BINNUM_PEAK] [--FC_thresh FC_THRESH] --output_dir OUTPUT_DIR --chrom_sizes CHROM_SIZES
                        [--chroms CHROMS [CHROMS ...]]
    
    Calculate local correlation and enrichment significance between two BedGraph tracks.
    
    options:  
    -h, --help                                          `Show this help message and exit`  
    --track1 TRACK1                                     `First input BedGraph file`  
    --track2 TRACK2                                     `Second input BedGraph file` 
    --output_dir OUTPUT_DIR                             `Output directory` 
    --chrom_sizes CHROM_SIZES                           `Path to the chromosome sizes file`
    --method {locP_and_ES,locP_and_ES}                  `Methods for calculate weighted local correlation and enrichment significance (default: locP_and_ES). P: pearson correlation; S: spearman correlation`  
    --FDR                                               `Use Benjamini–Hochberg FDR-corrected q-values instead of raw P-values` 
    --bin_size BIN_SIZE                                 `Size of each bin (default: 200)` 
    --binNum_window BIN_NUM_WINDOW                      `Number of bins in the sliding window (default: 11)`  
    --step STEP                                         `Step size for the sliding window (default: 1)`  
    --binNum_peak BIN_NUM_PEAK                          `Number of bins of the peak for ES (default: 11). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended`  
    --FC_thresh FC_THRESH                               `Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to the FC_thresh in the vocalno plot. Wald, a statistical value following a normal distribution here, euqal to log1.5(track1/track2) / SE can be used to calculate the p value, whose -log10 represents for ES here`
    --percentile PERCENTILE                             `Percentile for floor correction of low-coverage bins (default: 5). High percentile such as 90 or 95 is recommended, when tracks mainly contains some high sharp peaks, while small percentile like 5 is recommended when tracks mainly contain broad and relatively low peaks`   
    --chroms CHROMS [CHROMS ...]                        `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`
    
    Usage Example 1:
        localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --output_dir ./calc_out --chrom_sizes hg19.chrom.sizes --method locP_and_ES  --bin_size 200 --binNum_window 11 --step 1 --binNum_peak 11 --FC_thresh 1.5 --percentile 5 --chroms chr1 chr2
    
    Usage Example 2:
        localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --output_dir ./calc_out --chrom_sizes hg19.chrom.sizes --method locP_and_ES --FDR --bin_size 200 --binNum_window 11 --step 1 --binNum_peak 11 --FC_thresh 1.5 --percentile 5 --chroms all

### findreg
```bash
localfinder findreg -h
```
    usage: localfinder findreg [-h] --track_E TRACK_E --track_C TRACK_C --output_dir OUTPUT_DIR --chrom_sizes CHROM_SIZES [--p_thresh P_THRESH] [--min_regionSize MIN_REGIONSIZE] [--chroms CHROMS [CHROMS ...]]
    
    Identify genomic regions that show significant differences in correlation and enrichment.
    
    options:  
    -h, --help                                          `show this help message and exit`  
    --track_E TRACK_E                                   `Enrichment Significance BedGraph file`  
    --track_C TRACK_C                                   `Local Correlation BedGraph file`  
    --output_dir OUTPUT_DIR                             `Output directory for significant regions` 
    --chrom_sizes CHROM_SIZES                           `Path to the chromosome sizes file` 
    --p_thresh P_THRESH                                 `P-value threshold (default: 0.05)`  
    --binNum_thresh BINNUM_THRESH                       `Min consecutive significant bins per region (default: 2)`    
    --chroms CHROMS [CHROMS ...]                        `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`  

    
    Usage Example 1:
        localfinder findreg --track_E track_ES.bedgraph --track_C track_hmC.bedgraph --output_dir ./findreg_out --chrom_sizes hg19.chrom.sizes --p_thresh 0.05 --binNum_thresh 2 --chroms chr1 chr2
    
    Usage Example 2:
        localfinder findreg --track_E track_ES.bedgraph --track_C track_hmC.bedgraph --output_dir ./findreg_out --chrom_sizes hg19.chrom.sizes --p_thresh 0.05 --binNum_thresh 2 --chroms all

### pipeline
```bash
localfinder pipeline -h
```
    usage: localfinder pipeline [-h] --input_files INPUT_FILES [INPUT_FILES ...] [--output_dir OUTPUT_DIR] --chrom_sizes CHROM_SIZES [--bin_size BIN_SIZE]
                            [--method {locP_and_ES,locS_and_ES}] [--FDR] [--binNum_window BINNUM_WINDOW] [--binNum_peak BINNUM_PEAK] [--step STEP] [--percentile PERCENTILE]
                            [--FC_thresh FC_THRESH] [--p_thresh P_THRESH] [--binNum_thresh BINNUM_THRESH] [--chroms CHROMS [CHROMS ...]]
    
    Run all steps of the localfinder pipeline sequentially.
    
    options:  
    -h, --help                                          `Show this help message and exit`  
    --input_files INPUT_FILES [INPUT_FILES ...]         `Input BigWig files to process`  
    --output_dir OUTPUT_DIR                             `Output directory for all results`  
    --chrom_sizes CHROM_SIZES                           `Path to the chromosome sizes file`   
    --method {locP_and_ES,locS_and_ES}                  `Method for calculate local correlation and enrichment significance (default: locP_and_ES)` 
    --FDR                                               `Use Benjamini–Hochberg FDR-corrected q-values instead of raw P-values`
    --bin_size BIN_SIZE                                 `Size of each bin for binning tracks (default: 200bp)` 
    --binNum_window BIN_NUM_WINDOW                      `Number of bins in the sliding window (default: 11)`  
    --step STEP                                         `Step size for the sliding window (default: 1)`  
    --binNum_peak BIN_NUM_PEAK                          `Number of bins of the peak for ES (default: 11). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended`  
    --FC_thresh FC_THRESH                               `Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to the FC_thresh in the vocalno plot. Wald, a statistical value following a normal distribution here, euqal to log1.5(track1/track2) / SE can be used to calculate the p value, whose -log10 represents for ES here`
    --percentile PERCENTILE                             `Percentile for floor correction of low-coverage bins (default: 5). High percentile such as 90 or 95 is recommended, when tracks mainly contains some high sharp peaks, while small percentile like 5 is recommended when tracks mainly contain broad and relatively low peaks`
    --p_thresh P_THRESH                                 `P-value threshold (default: 0.05)`
    --binNum_thresh BINNUM_THRESH                       `Min consecutive significant bins per region (default: 2)` 
    --chroms CHROMS [CHROMS ...]                        `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`  

    
    Usage Example 1:
        localfinder pipeline --input_files track1.bedgraph track2.bedgraph --output_dir ./pipeline_out --chrom_sizes hg19.chrom.sizes --method locP_and_ES --bin_size 200 --binNum_window 11 --step 1 --binNum_peak 5 --FC_thresh 1.5 --percentile 5 --p_thresh 0.05 --binNum_thresh 2 --chroms chr1 chr2
    
    Usage Example 2:
        localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./pipeline_out --chrom_sizes hg19.chrom.sizes --method locP_and_ES --FDR --bin_size 200 --binNum_window 11 --step 1 --binNum_peak 5 --FC_thresh 1.5 --percentile 5 --p_thresh 0.05 --binNum_thresh 2 --chroms all

### viz
```bash
localfinder viz -h
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
git clone git@github.com:astudentfromsustech/localfinder.git
```

Run the examples under localfinder/tests/ (scripts have been preprared in tests folder)  
