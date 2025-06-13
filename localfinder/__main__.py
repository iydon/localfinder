# File: __main__.py

import argparse
import sys
import importlib.metadata
import argcomplete  # Import argcomplete for auto-completion
import textwrap  # Import textwrap for dedent
from argparse import RawDescriptionHelpFormatter  # Import RawDescriptionHelpFormatter

from localfinder.commands.bin import main as bin_tracks_main
from localfinder.commands.calc import main as calc_corr_main
from localfinder.commands.findreg import main as find_regions_main
from localfinder.commands.viz import main as visualize_main
from localfinder.pipeline import run_pipeline  # Import from pipeline.py

def main():
    # Retrieve package version
    try:
        version = importlib.metadata.version("localfinder")
    except importlib.metadata.PackageNotFoundError:
        version = "0.0.0"  # Fallback version

    # Create the top-level parser
    parser = argparse.ArgumentParser(
        prog='localfinder',
        description='localfinder: A tool calculating weighted local correlation and enrichment significance of two tracks and finding significantly different genomic regions. (GitHub: https://github.com/astudentfromsustech/localfinder)'
    )
    parser.add_argument('--version', '-V', action='version',
                        version=f'localfinder {version}',
                        help='Show program\'s version number and exit.')

    # Create subparsers for subcommands
    subparsers = parser.add_subparsers(dest='command', title='Subcommands',
                                       description='Valid subcommands',
                                       help='Additional help', metavar='')

    # Subcommand: bin (alias: bin_tracks)
    parser_bin = subparsers.add_parser(
        'bin',
        aliases=['bin_tracks'],
        help='Convert input files into bins with BedGraph format.',
        description='Bin genomic tracks into fixed-size bins and output BedGraph format.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder bin --input_files track1.bw track2.bw --output_dir ./binned_tracks --bin_size 200 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

            Usage Example 2:
                localfinder bin --input_files track1.bigwig track2.bigwig --output_dir ./binned_tracks --bin_size 200 --chrom_sizes hg19.chrom.sizes --chroms all
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )
    parser_bin.add_argument('--input_files', nargs='+', required=True,
                            help='Input files in BigWig/BedGraph/BAM/SAM format.')
    parser_bin.add_argument('--output_dir', required=True,
                            help='Output directory for binned data.')
    parser_bin.add_argument('--bin_size', type=int, default=200,
                            help='Size of each bin (default: 200).')
    parser_bin.add_argument('--chrom_sizes', type=str, required=True,
                            help='Path to the chromosome sizes file.')
    parser_bin.add_argument('--chroms', nargs='+', default=['all'],
                            help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".')
    parser_bin.set_defaults(func=bin_tracks_main)

    # Subcommand: calc (alias: calculate_localCorrelation_and_enrichmentSignificance)
    parser_calc = subparsers.add_parser(
        'calc',
        aliases=['calculate_localCorrelation_and_enrichmentSignificance'],
        help='Calculate weighted local correlation and enrichment significance between tracks.',
        description='Calculate weighted local correlation and enrichment significance between two BedGraph tracks.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --method locP_and_ES  --binNum_window 11 --step 1 --binNum_peak 11 --FC_thresh 1.5 --percentile 5 --output_dir ./results --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

            Usage Example 2:
                localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --method locP_and_ES  --binNum_window 11 --step 1 --binNum_peak 11 --FC_thresh 1.5 --percentile 5  --output_dir ./results --chrom_sizes hg19.chrom.sizes --chroms all
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )
    parser_calc.add_argument('--track1', required=True,
                             help='First input BedGraph file.')
    parser_calc.add_argument('--track2', required=True,
                             help='Second input BedGraph file.')
    parser_calc.add_argument('--method', choices=[
        'locP_and_ES',
        'locS_and_ES',
    ], default='locP_and_ES',  # Set default to one of the aliases
        help='Method for calculate weighted local correlation and enrichment significance (default: locP_and_ES)')
    parser_calc.add_argument('--binNum_window', type=int, default=11,
                             help='Number of bins in the sliding window (default: 11).')
    parser_calc.add_argument('--step', type=int, default=1,
                             help='Step size for the sliding window (default: 1).')
    parser_calc.add_argument('--percentile', type=int, default=5,
                             help='Percentile for floor correction of low-coverage bins (default: 5). High percentile such as 90 or 95 is recommended, when tracks mainly contains some high sharp peaks, while small percentile like 5 is recommended when tracks mainly contain broad and relatively low peaks')
    # NEW: small window size for the 3-bin ES calculation
    parser_calc.add_argument('--binNum_peak',type=int, default=11,
                            help='Number of bins of the peak for ES (default: 11). When the peak is around 400bp and the bin_size=200bp, binNum_peak=2 is recommended')
    # NEW: FC threshold base for the log-fold enrichment
    parser_calc.add_argument('--FC_thresh', type=float, default=1.5,
                            help='Fold-change threshold used as log base in enrichment (default: 1.5). When FC_thresh=1.5, the null hypothesis is that log1.5(track1/track2)=0, which is quite similar to the FC_thresh in the vocalno plot. Wald, a statistical value following a normal distribution here, euqal to log1.5(track1/track2) / SE can be used to calculate the p value, whose -log10 represents for ES here')
    parser_calc.add_argument('--output_dir', required=True,
                             help='Output directory for results.')
    parser_calc.add_argument('--chrom_sizes', type=str, required=True,
                             help='Path to the chromosome sizes file.')
    parser_calc.add_argument('--chroms', nargs='+', default=['all'],
                            help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".')
    parser_calc.set_defaults(func=calc_corr_main)

    # Subcommand: findreg (alias: find_significantly_different_regions)
    parser_find = subparsers.add_parser(
        "findreg", aliases=["find_significantly_different_regions"],
        help="Merge consecutive significant bins into regions. And find significantly different regions from ES & hmwC tracks.",
        description="Merge consecutive significant bins into regions. And find significantly different regions from ES & hmwC tracks.",
        epilog=textwrap.dedent("""\
            Example 1:
              localfinder findreg --track_E track_ES.bedgraph --track_C track_hmwC.bedgraph --output_dir ./findreg_out --p_thresh 0.05 --binNum_thresh 2 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2
            
            Example 2:
              localfinder findreg --track_E track_ES.bedgraph --track_C track_hmwC.bedgraph --output_dir ./findreg_out --p_thresh 0.05 --binNum_thresh 2 --chrom_sizes hg19.chrom.sizes --chroms all
        """),
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser_find.add_argument("--track_E", required=True, help="track_ES.bedgraph")
    parser_find.add_argument("--track_C", required=True, help="track_hmwC.bedgraph")
    parser_find.add_argument("--output_dir", required=True)
    # thresholds
    parser_find.add_argument("--p_thresh", type=float, default=0.05,    ### <<< NEW
                             help="P-value threshold (default 0.05)")
    parser_find.add_argument("--binNum_thresh", type=int, default=2,    ### <<< NEW
                             help="Min consecutive significant bins per region")
    parser_find.add_argument("--chroms", nargs="+", default=["all"])
    parser_find.add_argument("--chrom_sizes", required=True)
    parser_find.set_defaults(func=find_regions_main)

    # Subcommand: viz (alias: visualize_tracks_or_scatters)
    parser_visualize = subparsers.add_parser(
        'viz',
        aliases=['visualize_tracks'],
        help='Visualize genomic tracks.',
        description='Visualize genomic tracks.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder viz --input_files track1.bedgraph track2.bedgraph --output_file output.html --method plotly --region chr1 1000000 2000000 --colors blue red

            Usage Example 2:
                localfinder viz --input_files track1.bedgraph track2.bedgraph --output_file output.png --method pyGenomeTracks --region chr1 1000000 2000000 --colors
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )
    parser_visualize.add_argument('--input_files', nargs='+', required=True,
                                  help='Input BedGraph files to visualize.')
    parser_visualize.add_argument('--output_file', required=True,
                                  help='Output visualization file (e.g., PNG, HTML).')
    parser_visualize.add_argument('--method', choices=['pyGenomeTracks', 'plotly'], required=True,
                                  help='Visualization method to use.')
    parser_visualize.add_argument('--region', nargs=3, metavar=('CHROM', 'START', 'END'),
                                  help='Region to visualize in the format: CHROM START END (e.g., chr20 1000000 2000000).')
    parser_visualize.add_argument('--colors', nargs='+',
                                  help='Colors for the tracks (optional).')
    parser_visualize.set_defaults(func=visualize_main)

    # Subcommand: pipeline
    parser_pipe = subparsers.add_parser(
        'pipeline',
        help='Run the full pipeline.',
        description='Run all steps of the localfinder pipeline sequentially.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder pipeline --input_files track1.bedgraph track2.bedgraph --output_dir ../results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --binNum_window 11 --binNum_peak 11 --step 1 --percentile 5 --FC_thresh 1.5 --p_thresh 0.05 --binNum_thresh 2 --chroms chr1 chr2

            Usage Example 2:
                localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ../results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --binNum_window 11 --binNum_peak 11 --step 1 --percentile 5 --FC_thresh 1.5 --p_thresh 0.05 --binNum_thresh 2 --chroms all
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )

    parser_pipe.add_argument("--input_files", nargs="+", required=True)
    parser_pipe.add_argument("--output_dir", default="./output_pipeline")
    parser_pipe.add_argument("--chrom_sizes", required=True)
    parser_pipe.add_argument("--bin_size", type=int, default=200)
    # calc options forwarded
    parser_pipe.add_argument("--method", choices=["locP_and_ES", "locS_and_ES"],
                             default="locP_and_ES")
    parser_pipe.add_argument("--binNum_window", type=int, default=11)
    parser_pipe.add_argument("--binNum_peak",   type=int, default=11)
    parser_pipe.add_argument("--step", type=int, default=1)
    parser_pipe.add_argument("--percentile", type=int, default=5)
    parser_pipe.add_argument("--FC_thresh", type=float, default=1.5)
    # findreg thresholds forwarded
    parser_pipe.add_argument("--p_thresh", type=float, default=0.05,    ### <<< NEW
                             help="P-value threshold for ES significance")
    parser_pipe.add_argument("--binNum_thresh", type=int, default=2,    ### <<< NEW
                             help="Min consecutive significant bins per region")
    parser_pipe.add_argument('--chroms', nargs='+', default=['all'],
                            help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".')
    parser_pipe.set_defaults(func=run_pipeline)

    # Enable auto-completion
    argcomplete.autocomplete(parser)

    # Parse the arguments
    args = parser.parse_args()

    # Execute the appropriate function based on the subcommand
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == '__main__':
    main()
