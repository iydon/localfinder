# On Apple Silicon (osx-arm64) Bioconda doesn’t yet build the UCSC “bigWigToBedGraph” utility, so conda/mamba can’t find ucsc-bigwigtobedgraph
# Users need to download the tool, and export its PATH to the conda environment you work with

# download the bigWigToBedGraph
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigWigToBedGraph

# enter the conda env localfinder (or other name)
if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "❌  Please 'conda activate localfinder' (or your env) first."
  exit 1
fi

mkdir -p "$CONDA_PREFIX/bin"
cp ./bigWigToBedGraph $CONDA_PREFIX/bin/

# make is executable
chmod +x $CONDA_PREFIX/bin/bigWigToBedGraph

# chech whether bigWigToBedGraph can be used in the conda env
which bigWigToBedGraph # output should be $CONDA_PREFIX/bin/bigWigToBedGraph
