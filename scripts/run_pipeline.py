#!/usr/bin/env python3
"""
End-to-end analysis pipeline for viral linear complete genomes.

This script orchestrates the complete analysis workflow:
1. Load and organize metadata
2. Filter for exclusively linear vOTUs
3. Perform quality-based filtering
4. Analyze terminal sequences
5. Apply taxonomy-based filtering
6. Generate summary statistics

Usage:
    python scripts/run_pipeline.py --config config/config.yaml [--output-dir ./results]
"""

import argparse
import sys
from pathlib import Path
import yaml
import pandas as pd
import logging

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from data_loading import (
    load_metadata, load_representative_uvigs, load_excluded_sequences,
    extract_topology_quality_subsets, extract_quality_tiers, get_votu_representatives,
    save_results
)
from filtering import (
    filter_exclusive_linear_votus, remove_nested_votus,
    filter_by_length_difference, filter_by_group_size,
    filter_linear_uvigs_comprehensive, filter_by_min_genome_length,
    print_filter_summary
)
from taxonomy import (
    extract_ictv_and_host_class, add_genome_type_column,
    get_min_genome_lengths, get_class_summary
)
from terminal_analysis import (
    flag_representative_uvigs, compute_uvig_mismatches_per_member,
    filter_by_passing_members, summarize_passing_uvigs
)


def setup_logging(config: dict) -> logging.Logger:
    """Setup logging configuration."""
    log_level = getattr(logging, config.get('logging', {}).get('level', 'INFO'))
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    logging.basicConfig(level=log_level, format=log_format)
    logger = logging.getLogger(__name__)

    return logger


def main(config_path: str, output_dir: str = None):
    """
    Run complete analysis pipeline.

    Parameters
    ----------
    config_path : str
        Path to YAML configuration file.
    output_dir : str, optional
        Override output directory from config.
    """
    # Load configuration
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Setup logging
    logger = setup_logging(config)
    logger.info("=" * 70)
    logger.info("VIRAL LINEAR COMPLETE GENOMES ANALYSIS PIPELINE")
    logger.info("=" * 70)

    # Setup output directories
    if output_dir:
        config['output']['results_dir'] = output_dir

    for dir_key in ['results_dir', 'figures_dir', 'tables_dir', 'filtered_uvigs_dir']:
        Path(config['output'][dir_key]).mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory: {config['output']['results_dir']}")

    # ========================================================================
    # STEP 1: Load data
    # ========================================================================
    logger.info("\n[Step 1] Loading metadata...")
    df_filtered = load_metadata(config['data']['input_metadata'])
    logger.info(f"  Loaded {len(df_filtered)} total records")

    representative = load_representative_uvigs(config['data']['header_repr'])
    logger.info(f"  Found {len(representative)} representative UViGs")

    excluded = load_excluded_sequences(config['data']['exclude_nested'])
    logger.info(f"  Will exclude {len(excluded)} nested UViGs")

    # ========================================================================
    # STEP 2: Extract topology/quality subsets
    # ========================================================================
    logger.info("\n[Step 2] Extracting topology and quality subsets...")
    subsets = extract_topology_quality_subsets(df_filtered)

    for name, subset in subsets.items():
        n_votu = subset['votu'].nunique()
        logger.info(f"  {name}: {n_votu:,} vOTUs")

    # ========================================================================
    # STEP 3: Filter for exclusively linear vOTUs
    # ========================================================================
    logger.info("\n[Step 3] Identifying exclusively linear vOTUs...")
    cleaned_linear, n_exclusive = filter_exclusive_linear_votus(subsets['all_linear'])
    logger.info(f"  Found {n_exclusive:,} exclusively linear vOTUs")

    cleaned_linear, n_removed = remove_nested_votus(cleaned_linear, excluded)
    logger.info(f"  Removed {n_removed:,} vOTUs with nested fragments")

    # ========================================================================
    # STEP 4: Extract quality tiers
    # ========================================================================
    logger.info("\n[Step 4] Extracting quality tiers...")
    quality_tiers = extract_quality_tiers(cleaned_linear)

    for tier_name, tier_df in quality_tiers.items():
        n_votu = tier_df['votu'].nunique()
        logger.info(f"  {tier_name}: {n_votu:,} vOTUs")

    # ========================================================================
    # STEP 5: Perform comprehensive filtering for each quality tier
    # ========================================================================
    logger.info("\n[Step 5] Performing comprehensive filtering by quality tier...")

    correspondence = pd.merge(representative, df_filtered)
    max_bp = config['filtering']['max_bp_difference']
    min_uvig = config['filtering']['min_uvigs_per_votu']

    filtered_sets = {}
    for tier_name, tier_df in quality_tiers.items():
        logger.info(f"\n  Filtering {tier_name}...")
        filtered, stats = filter_linear_uvigs_comprehensive(
            cleaned_linear, tier_df, correspondence,
            max_bp_diff=max_bp,
            min_uvigs_per_votu=min_uvig
        )
        filtered_sets[tier_name] = filtered
        logger.info(f"    Final: {stats['after_group_size_filter']['votu']:,} vOTUs")

    # ========================================================================
    # STEP 6: Apply taxonomy filtering
    # ========================================================================
    logger.info("\n[Step 6] Applying taxonomy-based filtering...")

    # Combine all filtered sets
    combined_filtered = pd.concat(list(filtered_sets.values()), ignore_index=True)
    logger.info(f"  Combined total: {combined_filtered['votu'].nunique():,} vOTUs")

    # Extract taxonomy
    combined_filtered = extract_ictv_and_host_class(combined_filtered)
    combined_filtered = add_genome_type_column(combined_filtered)

    # Filter by minimum genome length
    min_lengths = get_min_genome_lengths()
    combined_filtered = filter_by_min_genome_length(combined_filtered, min_lengths)
    logger.info(f"  After length filter: {combined_filtered['votu'].nunique():,} vOTUs")

    # ========================================================================
    # STEP 7: Terminal sequence analysis
    # ========================================================================
    logger.info("\n[Step 7] Analyzing terminal sequences...")
    logger.info("  (Skipped in this example - load terminal data as needed)")

    # ========================================================================
    # STEP 8: Summary statistics
    # ========================================================================
    logger.info("\n[Step 8] Generating summary statistics...")

    class_summary = get_class_summary(combined_filtered)
    logger.info("\n  vOTU counts by viral class:")
    for virus_class, count in class_summary.head(10).items():
        logger.info(f"    {virus_class}: {count:,}")

    genome_type_summary = combined_filtered['Genome_type'].value_counts()
    logger.info("\n  vOTU counts by genome type:")
    for genome_type, count in genome_type_summary.items():
        logger.info(f"    {genome_type}: {count:,}")

    source_summary = combined_filtered['source'].value_counts()
    logger.info("\n  vOTU counts by source:")
    for source, count in source_summary.items():
        logger.info(f"    {source}: {count:,}")

    # ========================================================================
    # STEP 9: Save results
    # ========================================================================
    logger.info("\n[Step 9] Saving results...")

    save_results(
        combined_filtered,
        Path(config['output']['tables_dir']) / 'final_filtered_linear_complete.tsv'
    )

    # Save by quality tier
    for tier_name, tier_df in filtered_sets.items():
        save_results(
            tier_df,
            Path(config['output']['tables_dir']) / f'linear_complete_{tier_name}.tsv'
        )

    # Save summary statistics
    summary_stats = {
        'metric': ['Total vOTUs', 'Total UViGs', 'Unique Sources', 'Unique Classes'],
        'count': [
            combined_filtered['votu'].nunique(),
            len(combined_filtered),
            combined_filtered['source'].nunique(),
            combined_filtered['class'].nunique()
        ]
    }
    summary_df = pd.DataFrame(summary_stats)
    save_results(
        summary_df,
        Path(config['output']['tables_dir']) / 'summary_statistics.tsv'
    )

    logger.info("\n" + "=" * 70)
    logger.info("PIPELINE COMPLETED SUCCESSFULLY")
    logger.info("=" * 70)
    logger.info(f"Results saved to: {config['output']['results_dir']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Viral Linear Complete Genomes Analysis Pipeline'
    )
    parser.add_argument(
        '--config',
        required=True,
        help='Path to YAML configuration file'
    )
    parser.add_argument(
        '--output-dir',
        help='Override output directory from config'
    )

    args = parser.parse_args()

    try:
        main(args.config, args.output_dir)
    except Exception as e:
        logging.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)
