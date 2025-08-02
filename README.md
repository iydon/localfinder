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
