"""
Example: Basic filtering workflow

This script demonstrates the most common use case:
loading data, filtering by quality, and analyzing results.

Run from repository root:
    python examples/01_basic_filtering.py
"""

import sys
from pathlib import Path

# Setup paths
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'src'))

import pandas as pd
import yaml
from data_loading import (
    load_metadata, load_representative_uvigs, load_excluded_sequences,
    extract_topology_quality_subsets, extract_quality_tiers
)
from filtering import (
    filter_exclusive_linear_votus, remove_nested_votus,
    filter_linear_uvigs_comprehensive, filter_by_min_genome_length
)
from taxonomy import (
    extract_ictv_and_host_class, add_genome_type_column, get_min_genome_lengths
)


def main():
    """Run basic filtering example."""

    # Load configuration
    config_path = PROJECT_ROOT / 'config' / 'config.yaml'
    with open(config_path) as f:
        config = yaml.safe_load(f)

    print("=" * 70)
    print("BASIC FILTERING WORKFLOW")
    print("=" * 70)

    # ====================================================================
    # STEP 1: Load data
    # ====================================================================
    print("\n[Step 1] Loading metadata...")
    df_filtered = load_metadata(config['data']['input_metadata'])
    print(f"✓ Loaded {len(df_filtered):,} records")
    print(f"  Columns: {', '.join(df_filtered.columns[:5])}...")

    representative = load_representative_uvigs(config['data']['header_repr'])
    print(f"✓ Loaded {len(representative):,} representative UViGs")

    excluded = load_excluded_sequences(config['data']['exclude_nested'])
    print(f"✓ Will exclude {len(excluded):,} nested sequences")

    # ====================================================================
    # STEP 2: Extract subsets by topology
    # ====================================================================
    print("\n[Step 2] Extracting topology subsets...")
    subsets = extract_topology_quality_subsets(df_filtered)

    for name, subset in subsets.items():
        n_votu = subset['votu'].nunique()
        n_uvig = len(subset)
        print(f"  {name:20s}: {n_votu:7,} vOTUs ({n_uvig:,} UViGs)")

    # ====================================================================
    # STEP 3: Filter for exclusively linear topology
    # ====================================================================
    print("\n[Step 3] Identifying exclusively linear vOTUs...")
    cleaned_linear, n_exclusive = filter_exclusive_linear_votus(subsets['all_linear'])
    print(f"  Found {n_exclusive:,} exclusively linear vOTUs")

    cleaned_linear, n_removed = remove_nested_votus(cleaned_linear, excluded)
    print(f"  Removed {len(excluded):,} nested sequences")
    print(f"  Remaining: {cleaned_linear['votu'].nunique():,} vOTUs")

    # ====================================================================
    # STEP 4: Extract quality tiers
    # ====================================================================
    print("\n[Step 4] Extracting quality tiers...")
    quality_tiers = extract_quality_tiers(cleaned_linear)

    for tier_name, tier_df in quality_tiers.items():
        n_votu = tier_df['votu'].nunique()
        print(f"  {tier_name:25s}: {n_votu:7,} vOTUs")

    # ====================================================================
    # STEP 5: Filter each tier (example: high-quality only)
    # ====================================================================
    print("\n[Step 5] Filtering high-quality tier...")

    correspondence = pd.merge(representative, df_filtered)
    hq_tier = quality_tiers['high_quality']

    filtered_hq, stats = filter_linear_uvigs_comprehensive(
        cleaned_linear, hq_tier, correspondence,
        max_bp_diff=config['filtering']['max_bp_difference'],
        min_uvigs_per_votu=config['filtering']['min_uvigs_per_votu']
    )

    print("  Filtering steps:")
    for step, counts in stats.items():
        print(f"    {step:30s}: {counts['votu']:6,} vOTUs ({counts['uvig']:7,} UViGs)")

    # ====================================================================
    # STEP 6: Apply taxonomy filtering
    # ====================================================================
    print("\n[Step 6] Applying taxonomy filtering...")

    # Extract taxonomy information
    filtered_hq = extract_ictv_and_host_class(filtered_hq)
    filtered_hq = add_genome_type_column(filtered_hq)

    print(f"  Before length filter: {filtered_hq['votu'].nunique():,} vOTUs")

    # Filter by minimum genome length
    min_lengths = get_min_genome_lengths()
    filtered_hq = filter_by_min_genome_length(filtered_hq, min_lengths)

    print(f"  After length filter: {filtered_hq['votu'].nunique():,} vOTUs")

    # ====================================================================
    # STEP 7: Summary statistics
    # ====================================================================
    print("\n[Step 7] Summary statistics")

    print("\n  Viral classes (top 10):")
    class_counts = filtered_hq['class'].value_counts().head(10)
    for virus_class, count in class_counts.items():
        print(f"    {virus_class:30s}: {count:5,} UViGs")

    print("\n  Genome types:")
    genome_counts = filtered_hq['Genome_type'].value_counts()
    for gtype, count in genome_counts.items():
        print(f"    {gtype:30s}: {count:5,} UViGs")

    print("\n  Data sources:")
    source_counts = filtered_hq['source'].value_counts()
    for source, count in source_counts.items():
        print(f"    {source:30s}: {count:5,} UViGs")

    # ====================================================================
    # STEP 8: Save results
    # ====================================================================
    print("\n[Step 8] Saving results...")

    output_dir = PROJECT_ROOT / 'results'
    output_dir.mkdir(exist_ok=True)

    output_file = output_dir / 'example_filtered_high_quality.tsv'
    filtered_hq.to_csv(output_file, sep='\t', index=False)
    print(f"✓ Saved to {output_file}")

    # ====================================================================
    # Detailed statistics
    # ====================================================================
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nFinal dataset:")
    print(f"  Total UViGs: {len(filtered_hq):,}")
    print(f"  Unique vOTUs: {filtered_hq['votu'].nunique():,}")
    print(f"  Unique sources: {filtered_hq['source'].nunique()}")
    print(f"  Length range: {filtered_hq['length'].min():,} - {filtered_hq['length'].max():,} bp")
    print(f"  Mean length: {filtered_hq['length'].mean():,.0f} bp")
    print(f"  Median length: {filtered_hq['length'].median():,.0f} bp")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
