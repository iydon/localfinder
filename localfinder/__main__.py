# File: __main__.py

import argparse
import json
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
        description='localfinder: A tool calculating local correlation and enrichment significance of two tracks and finding significantly different genomic regions. (GitHub: https://github.com/astudentfromsustech/localfinder)'
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
        help='Calculate local correlation and enrichment significance between tracks.',
        description='Calculate local correlation and enrichment significance between two BedGraph tracks.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --method locP_and_ES --method_params '{"percentile": 5}' --bin_num 11 --step 1 --output_dir ./results --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

            Usage Example 2:
                localfinder calc --track1 track1.bedgraph --track2 track2.bedgraph --method locP_and_ES --method_params '{"percentile": 5}' --bin_num 11 --step 1 --output_dir ./results --chrom_sizes hg19.chrom.sizes --chroms all
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )
    parser_calc.add_argument('--track1', required=True,
                             help='First input BedGraph file.')
    parser_calc.add_argument('--track2', required=True,
                             help='Second input BedGraph file.')
    parser_calc.add_argument('--method', choices=[
        'locP_and_ES',
        'locWP_and_ES',
        'locS_and_ES',
        'locWS_and_ES',
        'locMI_and_ES'
    ], default='locP_and_ES',  # Set default to one of the aliases
        help='Method for calculate local correlation and enrichment significance (default: locP_and_ES)')
    parser_calc.add_argument('--method_params', type=str, default='{"percentile": 5}',
                             help='Parameters for the method in JSON format (default: {"percentile": 5}).')
    parser_calc.add_argument('--bin_num', type=int, default=11,
                             help='Number of bins in the sliding window (default: 11).')
    parser_calc.add_argument('--step', type=int, default=1,
                             help='Step size for the sliding window (default: 1).')
    parser_calc.add_argument('--output_dir', required=True,
                             help='Output directory for results.')
    parser_calc.add_argument('--chrom_sizes', type=str, required=True,
                             help='Path to the chromosome sizes file.')
    parser_calc.add_argument('--chroms', nargs='+', default=['all'],
                            help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".')
    parser_calc.set_defaults(func=calc_corr_main)

    # Subcommand: findreg (alias: find_significantly_different_regions)
    parser_find = subparsers.add_parser(
        'findreg',
        aliases=['find_significantly_different_regions'],
        help='Find significantly different regions between tracks.',
        description='Identify genomic regions that show significant differences in correlation and enrichment.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder findreg --track_E track_E.bedgraph --track_C track_C.bedgraph --output_dir ./results --min_regionSize 5 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chrom_sizes hg19.chrom.sizes --chroms chr1 chr2

            Usage Example 2:
                localfinder findreg --track_E track_E.bedgraph --track_C track_C.bedgraph --output_dir ./results --min_regionSize 5 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chrom_sizes hg19.chrom.sizes --chroms all
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )
    parser_find.add_argument('--track_E', required=True,
                             help='Enrichment Significance BedGraph file.')
    parser_find.add_argument('--track_C', required=True,
                             help='Local Correlation BedGraph file.')
    parser_find.add_argument('--output_dir', required=True,
                             help='Output directory for significant regions.')
    parser_find.add_argument('--min_regionSize', type=int, default=5,
                             help='Minimum number of consecutive bins to define a region (default: 5).')
    parser_find.add_argument('--E_upPercentile', type=float, default=75,
                             help='High percentile for enrichment (default: 75).')
    parser_find.add_argument('--E_lowPercentile', type=float, default=25,
                             help='Low percentile for enrichment (default: 25).')
    parser_find.add_argument('--C_upPercentile', type=float, default=75,
                             help='High percentile for correlation (default: 75).')
    parser_find.add_argument('--C_lowPercentile', type=float, default=25,
                             help='Low percentile for correlation (default: 25).')
    parser_find.add_argument('--chroms', nargs='+', default=['all'],
                            help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".')
    parser_find.add_argument('--chrom_sizes', type=str, required=True,
                             help='Path to the chromosome sizes file.')
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
    parser_pipeline = subparsers.add_parser(
        'pipeline',
        help='Run the full pipeline.',
        description='Run all steps of the localfinder pipeline sequentially.',
        epilog=textwrap.dedent('''\
            Usage Example 1:
                localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --bin_num 11 --step 1 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chroms chr1 chr2

            Usage Example 2:
                localfinder pipeline --input_files track1.bigwig track2.bigwig --output_dir ./results --chrom_sizes hg19.chrom.sizes --bin_size 200 --method locP_and_ES --bin_num 11 --step 1 --E_upPercentile 75 --E_lowPercentile 25 --C_upPercentile 75 --C_lowPercentile 25 --chroms all
            '''),
        formatter_class=RawDescriptionHelpFormatter  # Preserve formatting
    )
    # Define necessary arguments for pipeline
    parser_pipeline.add_argument('--input_files', nargs='+', required=True,
                                 help='Input BigWig files to process.')
    parser_pipeline.add_argument('--output_dir', type=str, default='./output_pipeline',
                                 help='Output directory for all results (default: ./output_pipeline).')
    parser_pipeline.add_argument('--chrom_sizes', type=str, required=True,
                                 help='Path to the chromosome sizes file.')
    parser_pipeline.add_argument('--bin_size', type=int, default=200,
                                 help='Size of each bin for binning tracks (default: 200bp)')
    parser_pipeline.add_argument('--method', choices=[
        'locP_and_ES',
        'locWP_and_ES',
        'locS_and_ES',
        'locWS_and_ES',
        'locMI_and_ES'
    ], default='locP_and_ES',
        help='Method for calculate local correlation and enrichment significance (default: locP_and_ES)')
    parser_pipeline.add_argument('--method_params', type=json.loads, default={"percentile": 5},
                                 help='Method-specific parameters in JSON format')
    parser_pipeline.add_argument('--bin_num', type=int, default=11,
                                 help='Number of bins in the sliding window (default: 11)')
    parser_pipeline.add_argument('--step', type=int, default=1,
                                 help='Step size for sliding window (default: 1)')
    parser_pipeline.add_argument('--E_upPercentile', type=float, default=75,
                                 help='Percentile threshold for high enrichment (default: 75)')
    parser_pipeline.add_argument('--E_lowPercentile', type=float, default=25,
                                 help='Percentile threshold for low enrichment (default: 25)')
    parser_pipeline.add_argument('--C_upPercentile', type=float, default=75,
                                 help='Percentile threshold for high correlation (default: 75)')
    parser_pipeline.add_argument('--C_lowPercentile', type=float, default=25,
                                 help='Percentile threshold for low correlation (default: 25)')
    parser_pipeline.add_argument('--chroms', nargs='+', default=['all'],
                                 help='Chromosomes to process (e.g., chr1 chr2). Defaults to "all".')
    parser_pipeline.add_argument('--min_regionSize', type=int, default=5,
                                 help='Minimum number of consecutive bins to define a region (default: 5).')
    parser_pipeline.set_defaults(func=run_pipeline)

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
