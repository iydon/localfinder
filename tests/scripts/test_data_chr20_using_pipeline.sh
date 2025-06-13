#!/bin/bash
export PYTHONUNBUFFERED=1

chrom_size="../annotations/hg19.chrom.sizes"
data_dir="../test_data_chr20"
result_dir="../result_test_data_chr20_using_pipeline"

method="locP_and_ES"


track_1="${data_dir}/track1_chr20.bedgraph"
track_2="${data_dir}/track2_chr20.bedgraph"


localfinder pipeline \
    --input_files ${track_1} ${track_2} \
    --output_dir  ${result_dir} \
    --chrom_sizes $chrom_size \
    --bin_size    200 \
    --method      locP_and_ES \
    --binNum_window 11 --binNum_peak 11 --step 1 \
    --percentile 5 --FC_thresh 1.5 \
    --p_thresh 0.05 --binNum_thresh 2 \
    --chroms chr20