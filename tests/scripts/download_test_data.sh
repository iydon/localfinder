#!/bin/bash

# Directory to save the downloaded files
OUTPUT_DIR="../data"

# URLs of the files to download
URLS=(
    "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E100-H3K4me1.pval.signal.bigwig"
    "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E071-H3K4me1.pval.signal.bigwig"
)

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Download each file into the directory
for URL in "${URLS[@]}"; do
    wget -P "$OUTPUT_DIR" "$URL"
done

echo "Download complete."

