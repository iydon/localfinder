#!/bin/bash
export PYTHONUNBUFFERED=1

# If there is no zstd tool, install it using conda install -c conda-forge zstd
zstdcat ../test_data_chr20/track1_chr20.bedgraph.zst > ../test_data_chr20/track1_chr20.bedgraph
zstdcat ../test_data_chr20/track2_chr20.bedgraph.zst > ../test_data_chr20/track2_chr20.bedgraph

chrom_size="../annotations/hg19.chrom.sizes"
data_dir="../test_data_chr20"
result_dir="../result_test_data_chr20_step_by_step"

method="locP_and_ES"
bin_size=200
percentile=5
p_thresh=0.05
binNum_thresh=2

track_1="${data_dir}/track1_chr20.bedgraph"
track_2="${data_dir}/track2_chr20.bedgraph"

bin_dir="${result_dir}/bin_tracks"
mkdir -p ${bin_dir}
localfinder bin --input_files ${track_1} ${track_2} --output_dir ${bin_dir} --bin_size ${bin_size} --chrom_sizes ${chrom_size} --chroms chr20



correlation_enrichment_dir="${result_dir}/correlation_enrichment"
mkdir -p ${correlation_enrichment_dir}

base1=$(basename "$track_1" .bedgraph)
bin_track_1="${result_dir}/bin_tracks/${base1}.binSize${bin_size}.bedgraph"
# echo ${bin_track_1}
base2=$(basename "$track_2" .bedgraph)
bin_track_2="${result_dir}/bin_tracks/${base2}.binSize${bin_size}.bedgraph"
# echo ${bin_track_2}

localfinder calc --track1 ${bin_track_1} --track2 ${bin_track_2} --output_dir ${correlation_enrichment_dir} \
--method ${method} --percentile ${percentile} --binNum_window 11 --binNum_peak 11 --step 1 --chrom_sizes ${chrom_size} --chroms chr20

significant_regions_dir="${result_dir}/significant_regions"
mkdir -p ${significant_regions_dir}

localfinder findreg --track_E ${correlation_enrichment_dir}/track_ES.bedgraph --track_C ${correlation_enrichment_dir}/track_hmC.bedgraph --output_dir ${significant_regions_dir} --p_thresh ${p_thresh} --binNum_thresh ${binNum_thresh} --chrom_size ${chrom_size} --chroms chr20