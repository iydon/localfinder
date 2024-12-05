# localfinder

A tool calculating local correlation and enrichment significance of two tracks and finding significantly different genomic regions.

## Installation Requirements

Before installing and using `localfinder`, please ensure that the following external tools are installed on your system:

- **bedtools**: Used for genomic interval operations.
  - Installation: [https://bedtools.readthedocs.io/en/latest/content/installation.html](https://bedtools.readthedocs.io/en/latest/content/installation.html)
  - conda install -c bioconda -c conda-forge bedtools 
  - mamba install -c bioconda -c conda-forge bedtools
- **ucsc-bigwigtobedgraph**: Used for converting BigWig files to BedGraph format.
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
usage: localfinder bin [-h] --input_files INPUT_FILES [INPUT_FILES ...] --output_dir OUTPUT_DIR [--bin_size BIN_SIZE] --chrom_sizes CHROM_SIZES [--chroms CHROMS [CHROMS ...]]

Bin genomic tracks into fixed-size bins and output BedGraph format.

options:  
  -h, --help `Show this help message and exit`  
  --input_files INPUT_FILES [INPUT_FILES ...] `Input files in BigWig/BedGraph/BAM/SAM format`    
  --output_dir OUTPUT_DIR `Output directory for binned data`  
  --bin_size BIN_SIZE `Size of each bin (default: 200)`  
  --chrom_sizes CHROM_SIZES `Path to the chromosome sizes file`  
  --chroms CHROMS [CHROMS ...] `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`

Usage Example 1:
    localfinder bin --input_files track1.bw track2.bw --output_dir ./binned_tracks --bin_size 200 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

Usage Example 2:
    localfinder bin --input_files track1.bigwig track2.bigwig --output_dir ./binned_tracks --bin_size 200 --chrom_sizes hg19.chrom.sizes --chroms all


### calc
```bash
localfinder calc -h
```
usage: localfinder calc [-h] --track1 TRACK1 --track2 TRACK2 [--method {locP_and_ES,locWP_and_ES,locS_and_ES,locWS_and_ES,locMI_and_ES}] [--method_params METHOD_PARAMS] [--bin_num BIN_NUM] [--step STEP] --output_dir OUTPUT_DIR --chrom_sizes CHROM_SIZES [--chroms CHROMS [CHROMS ...]]

Calculate local correlation and enrichment significance between two BedGraph tracks.

options:  
-h, --help            `Show this help message and exit`  
--track1 TRACK1       `First input BedGraph file`  
--track2 TRACK2       `Second input BedGraph file`  
--method {locP_and_ES,locWP_and_ES,locS_and_ES,locWS_and_ES,locMI_and_ES} `Methods for calculate local correlation and enrichment significance (default: locP_and_ES)`  
--method_params METHOD_PARAMS `Parameters for the method in JSON format (default: {"percentile": 5})`  
--bin_num BIN_NUM `Number of bins in the sliding window (default: 11)`  
--step STEP  `Step size for the sliding window (default: 1)`  
--output_dir OUTPUT_DIR `Output directory for results`  
--chrom_sizes CHROM_SIZES `Path to the chromosome sizes file`  
--chroms CHROMS [CHROMS ...] `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`

Usage Example 1:
    localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --method locP_and_ES --method_params '{"percentile": 5}' --bin_num 11 --step 1 --output_dir ./results --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

Usage Example 2:
    localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --method locP_and_ES --method_params '{"percentile": 5}' --bin_num 11 --step 1 --output_dir ./results --chrom_sizes hg19.chrom.sizes --chroms all  

### findreg
```bash
localfinder findreg -h
```
usage: localfinder findreg [-h] --track_E TRACK_E --track_C TRACK_C --output_dir OUTPUT_DIR [--min_regionSize MIN_REGIONSIZE] [--E_upPercentile E_UPPERCENTILE] [--E_lowPercentile E_LOWPERCENTILE] [--C_upPercentile C_UPPERCENTILE] [--C_lowPercentile C_LOWPERCENTILE] [--chroms CHROMS [CHROMS ...]]
                           --chrom_sizes CHROM_SIZES

Identify genomic regions that show significant differences in correlation and enrichment.

options:  
-h, --help            `show this help message and exit`  
--track_E TRACK_E     `Enrichment Significance BedGraph file`  
--track_C TRACK_C     `Local Correlation BedGraph file`  
--output_dir OUTPUT_DIR `Output directory for significant regions`  
--min_regionSize MIN_REGIONSIZE `Minimum number of consecutive bins to define a region (default: 5)`  
--E_upPercentile E_UPPERCENTILE `High percentile for enrichment (default: 75)`  
--E_lowPercentile E_LOWPERCENTILE `Low percentile for enrichment (default: 25)`  
--C_upPercentile C_UPPERCENTILE  `High percentile for correlation (default: 75)`  
--C_lowPercentile C_LOWPERCENTILE  `Low percentile for correlation (default: 25)`  
--chroms CHROMS [CHROMS ...] `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`  
--chrom_sizes CHROM_SIZES `Path to the chromosome sizes file`

Usage Example 1:
    localfinder findreg --track_E track_E.bedgraph --track_C track_C.bedgraph --output_dir ./results --min_regionSize 5 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

Usage Example 2:
    localfinder findreg --track_E track_E.bedgraph --track_C track_C.bedgraph --output_dir ./results --min_regionSize 5 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chrom_sizes hg19.chrom.sizes --chroms all

### pipeline
```bash
localfinder pipeline -h
```
usage: localfinder pipeline [-h] --input_files INPUT_FILES [INPUT_FILES ...] [--output_dir OUTPUT_DIR] --chrom_sizes CHROM_SIZES [--bin_size BIN_SIZE] [--method {locP_and_ES,locWP_and_ES,locS_and_ES,locWS_and_ES,locMI_and_ES}] [--method_params METHOD_PARAMS] [--bin_num BIN_NUM] [--step STEP]
                            [--E_upPercentile E_UPPERCENTILE] [--E_lowPercentile E_LOWPERCENTILE] [--C_upPercentile C_UPPERCENTILE] [--C_lowPercentile C_LOWPERCENTILE] [--chroms CHROMS [CHROMS ...]] [--min_regionSize MIN_REGIONSIZE]

Run all steps of the localfinder pipeline sequentially.

options:  
-h, --help            `Show this help message and exit`  
--input_files INPUT_FILES [INPUT_FILES ...] `Input BigWig files to process`  
--output_dir OUTPUT_DIR  `Output directory for all results (default: ./output_pipeline)`  
--chrom_sizes CHROM_SIZES  `Path to the chromosome sizes file`  
--bin_size BIN_SIZE   `Size of each bin for binning tracks (default: 200bp)`  
--method {locP_and_ES,locWP_and_ES,locS_and_ES,locWS_and_ES,locMI_and_ES} `Method for calculate local correlation and enrichment significance (default: locP_and_ES)`  
--method_params METHOD_PARAMS `Method-specific parameters in JSON format`  
--bin_num BIN_NUM     `Number of bins in the sliding window (default: 11)`  
--step STEP           `Step size for sliding window (default: 1)`  
--E_upPercentile E_UPPERCENTILE `Percentile threshold for high enrichment (default: 75)`  
--E_lowPercentile E_LOWPERCENTILE `Percentile threshold for low enrichment (default: 25)`  
--C_upPercentile C_UPPERCENTILE  `Percentile threshold for high correlation (default: 75)`  
--C_lowPercentile C_LOWPERCENTILE `Percentile threshold for low correlation (default: 25)`  
--chroms CHROMS [CHROMS ...] `Chromosomes to process (e.g., chr1 chr2). Defaults to "all"`  
--min_regionSize MIN_REGIONSIZE `Minimum number of consecutive bins to define a region (default: 5)`

Usage Example 1:
    localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --bin_num 11 --step 1 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chroms chr1 chr2

Usage Example 2:
    localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --bin_num 11 --step 1 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chroms all

### viz
```bash
localfinder viz -h
```
usage: localfinder viz [-h] --input_files INPUT_FILES [INPUT_FILES ...] --output_file OUTPUT_FILE --method {pyGenomeTracks,plotly} [--region CHROM START END] [--colors COLORS [COLORS ...]]

Visualize genomic tracks.

options:  
-h, --help            `show this help message and exit`  
--input_files INPUT_FILES [INPUT_FILES ...] `Input BedGraph files to visualize`  
--output_file OUTPUT_FILE  `Output visualization file (e.g., PNG, HTML)`  
--method {pyGenomeTracks,plotly} `Visualization method to use`  
--region CHROM START END `Region to visualize in the format: CHROM START END (e.g., chr20 1000000 2000000)`  
--colors COLORS [COLORS ...] `Colors for the tracks (optional)`

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
conda install -c bioconda -c conda-forge samtools bedtools ucsc-bigwigtobedgraph
pip install LocalFinder
```

Run the localfinder/tests/scripts/download_test_data.sh to download the test data
```bash
cd ./localfinder/tests/scripts
bash .download_test_data.sh
```

Next, we will run all the subcommands on the test data   
1. bin (bin_tracks) 
```bash
localfinder bin -h 
localfinder bin --input_files ./tests/data/E071-H3K4me1.pval.signal.bigwig ./tests/data/E100-H3K4me1.pval.signal.bigwig --output_dir ./tests/binned_tracks --chrom_sizes ./tests/annotations/hg19.chrom.sizes --chroms chr20 chr19
```

2. localfinder calc (calc_locCor_and_ES) 
```bash
localfinder calc -h 
localfinder calc --track1 ./tests/binned_tracks/E071-H3K4me1.pval.signal.binSize200.bedgraph --track2 ./tests/binned_tracks/E100-H3K4me1.pval.signal.binSize200.bedgraph --method locP_and_ES --output_dir ./tests/calc --chrom_sizes ./tests/annotations/hg19.chrom.sizes --chroms chr20 chr19
```

3. localfinder findreg (find_SDR)
```bash
localfinder findreg -h
localfinder findreg --track_E ./tests/calc/track_ES.bedgraph --track_C ./tests/calc/track_locCor.bedgraph --output_dir ./tests/findreg --chrom_sizes ./tests/annotations/hg19.chrom.sizes --chroms chr19 chr20
```

4. localfinder pipeline
```bash
localfinder pipeline -h 
localfinder pipeline --input_files ./tests/data/E071-H3K4me1.pval.signal.bigwig ./tests/data/E100-H3K4me1.pval.signal.bigwig --output_dir ./tests/pipeline --chrom_sizes ./tests/annotations/hg19.chrom.sizes --chroms chr20 chr19
```

5. localfinder viz (visualize_tracks)
```bash
localfinder viz -h
localfinder viz --input_files ./tests/binned_tracks/E071-H3K4me1.pval.signal.binSize200.bedgraph ./tests/binned_tracks/E100-H3K4me1.pval.signal.binSize200.bedgraph ./tests/calc/track_locCor.bedgraph ./tests/calc/track_ES.bedgraph --output_file ./tests/viz/test.html --method plotly --region chr20 1000000 2000000
```