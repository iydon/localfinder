#!/bin/bash
export PYTHONUNBUFFERED=1

# If there is no zstd tool, install it using conda install -c conda-forge zstd
zstdcat ../test_ChIP-seq_data_from_NG2011/m3134-pDEX-ChIP_allCombined.100.bw.zst > ../test_ChIP-seq_data_from_NG2011/m3134-pDEX-ChIP_allCombined.100.bw
zstdcat ../test_ChIP-seq_data_from_NG2011/m3134-DEX-ChIP_allCombined.100.bw.zst > ../test_ChIP-seq_data_from_NG2011/m3134-DEX-ChIP_allCombined.100.bw

chrom_size="../annotations/mm10.chrom.sizes"
data_dir="../test_ChIP-seq_data_from_NG2011"
result_dir="../result_ChIP-seq_data_from_NG2011_using_pipeline"

method="locP_and_ES"


track_1="${data_dir}/m3134-pDEX-ChIP_allCombined.100.bw"
track_2="${data_dir}/m3134-DEX-ChIP_allCombined.100.bw"


# After visulaization of the two tracks, we can find that most peaks are sharp peaks and most bins are with low values, so we choose large percentile 90 and binNum_peak as 5 (or 3)
localfinder pipeline \
    --input_files ${track_1} ${track_2} \
    --output_dir  ${result_dir} \
    --chrom_sizes $chrom_size \
    --bin_size    200 \
    --method      locP_and_ES \
    --binNum_window 11 --binNum_peak 5 --step 1 \
    --percentile 90 --FC_thresh 1.5 \
    --p_thresh 0.05 --binNum_thresh 2 \
    --chroms chr1 chr2
