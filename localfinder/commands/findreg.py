# File: commands/findreg.py

import os
from localfinder.utils import (                       # util helpers
    find_significantly_different_regions as util_find,  ### <<< CHANGED
    get_chromosomes_from_chrom_sizes
)

# ----------------------------------------------------------------------
# CLI entry-point (wired in __main__.py)
# ----------------------------------------------------------------------
def main(args):                                        ### <<< CHANGED
    # • pull CLI arguments
    track_ES_file   = args.track_E                     ### <<< CHANGED
    track_hmwC_file = args.track_C                     ### <<< CHANGED
    output_dir      = args.output_dir
    p_thresh        = args.p_thresh                    ### <<< NEW
    binNum_thresh   = args.binNum_thresh               ### <<< NEW
    chroms          = args.chroms
    chrom_sizes     = args.chrom_sizes

    os.makedirs(output_dir, exist_ok=True)

    # • expand “all” → actual chromosome list
    if chroms == ['all'] or chroms is None:
        chroms = get_chromosomes_from_chrom_sizes(chrom_sizes)
        print(f"'chroms' set to all chromosomes: {chroms}")
    else:
        print(f"'chroms' set to specified chromosomes: {chroms}")

    # • call the utility function (aliased as util_find)
    util_find(
        track_ES_file   = track_ES_file,               ### <<< CHANGED
        track_hmwC_file = track_hmwC_file,             ### <<< CHANGED
        output_dir      = output_dir,
        p_thresh        = p_thresh,                    ### <<< NEW
        binNum_thresh   = binNum_thresh,               ### <<< NEW
        chroms          = chroms,
        chrom_sizes     = chrom_sizes
    )