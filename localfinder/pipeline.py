# File: pipeline.py

import argparse, os, sys
from localfinder.commands.bin import main as bin_tracks_main
from localfinder.commands.calc import main as calc_corr_main
from localfinder.commands.findreg import main as find_regions_main  
from localfinder.utils import check_external_tools, get_chromosomes_from_chrom_sizes

def run_pipeline(args):
    # Ensure required external tools are available
    check_external_tools()

    # Step 1: Bin the tracks
    bin_output_dir = os.path.join(args.output_dir, 'binned_tracks')
    bin_args = argparse.Namespace(
        input_files=args.input_files,
        output_dir=bin_output_dir,
        bin_size=args.bin_size,
        chrom_sizes=args.chrom_sizes,
        chroms=args.chroms
    )
    bin_tracks_main(bin_args)

    # Step 2: Calculate correlation and enrichment
    binned_files = [
        os.path.join(bin_output_dir, os.path.basename(f).replace('.bam', f'.binSize{args.bin_size}.bedgraph')
                                   .replace('.sam', f'.binSize{args.bin_size}.bedgraph')
                                   .replace('.bedgraph', f'.binSize{args.bin_size}.bedgraph')
                                   .replace('.bigwig', f'.binSize{args.bin_size}.bedgraph')
                                   .replace('.bw', f'.binSize{args.bin_size}.bedgraph'))
        for f in args.input_files
    ]

    # Ensure that exactly two binned files are present
    if len(binned_files) < 2:
        print("Error: At least two binned files are required for correlation and enrichment calculation.")
        sys.exit(1)

    # **Modification Start**
    # If chroms is 'all', retrieve all chromosomes from chrom_sizes
    if args.chroms == ['all'] or args.chroms is None:
        chroms = get_chromosomes_from_chrom_sizes(args.chrom_sizes)
        print(f"'chroms' set to all chromosomes from chrom_sizes: {chroms}")
    else:
        chroms = args.chroms
        print(f"'chroms' set to specified chromosomes: {chroms}")
    # **Modification End**

    calc_output_dir = os.path.join(args.output_dir, 'correlation_enrichment')
    calc_args = argparse.Namespace(
        track1=binned_files[0],
        track2=binned_files[1],
        output_dir=calc_output_dir,
        method=args.method,
        percentile=args.percentile,
        binNum_window=args.binNum_window,  ### <<< CHANGED
        step=args.step,
        binNum_peak=args.binNum_peak,      ### <<< CHANGED
        FC_thresh=args.FC_thresh,          ### <<< CHANGED
        chroms=chroms,
        chrom_sizes=args.chrom_sizes  # **Pass chrom_sizes to calc.py**
    )
    calc_corr_main(calc_args)

    # ── 3. find significantly different regions ───────────────────────
    find_output_dir = os.path.join(args.output_dir, "significant_regions")
    find_args = argparse.Namespace(
        track_E   = os.path.join(calc_output_dir, "track_ES.bedgraph"),   ### FIX
        track_C   = os.path.join(calc_output_dir, "track_hmwC.bedgraph"), ### FIX
        output_dir      = find_output_dir,
        p_thresh        = args.p_thresh,        ### FIX
        binNum_thresh   = args.binNum_thresh,   ### FIX
        chroms          = chroms,
        chrom_sizes     = args.chrom_sizes,
    )
    find_regions_main(find_args)

    print("Pipeline completed successfully.")

