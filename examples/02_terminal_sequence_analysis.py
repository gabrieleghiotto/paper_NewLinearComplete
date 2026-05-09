"""
Example: Terminal sequence analysis

This script demonstrates terminal sequence conservation analysis:
- Load conserved start/end sequences
- Compare against representative sequences
- Identify orthologs based on sequence similarity

Run from repository root:
    python examples/02_terminal_sequence_analysis.py
"""

import sys
from pathlib import Path

# Setup paths
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / 'src'))

import pandas as pd
import yaml
from data_loading import load_metadata, load_representative_uvigs
from filtering import filter_by_min_genome_length
from taxonomy import (
    extract_ictv_and_host_class, add_genome_type_column, get_min_genome_lengths
)
from terminal_analysis import (
    flag_representative_uvigs, compute_uvig_mismatches_per_member,
    filter_by_passing_members, summarize_passing_uvigs
)


def main():
    """Run terminal sequence analysis example."""

    # Load configuration
    config_path = PROJECT_ROOT / 'config' / 'config.yaml'
    with open(config_path) as f:
        config = yaml.safe_load(f)

    print("=" * 70)
    print("TERMINAL SEQUENCE ANALYSIS")
    print("=" * 70)

    # ====================================================================
    # Load filtered dataset (would normally be from previous step)
    # ====================================================================
    print("\n[Step 1] Loading data...")

    # For this example, we'll use the saved result from basic filtering
    # In practice, you'd use your own filtered dataset
    filtered_file = PROJECT_ROOT / 'results' / 'example_filtered_high_quality.tsv'

    if not filtered_file.exists():
        print("⚠ Running basic filtering first...")
        from examples import example_01_basic_filtering
        example_01_basic_filtering.main()

    df_uvigs = pd.read_csv(filtered_file, sep='\t')
    print(f"✓ Loaded {len(df_uvigs):,} UViGs from {df_uvigs['votu'].nunique():,} vOTUs")

    # Load representative mapping
    representative = load_representative_uvigs(config['data']['header_repr'])
    print(f"✓ Loaded {len(representative):,} representatives")

    # ====================================================================
    # Load terminal sequences
    # ====================================================================
    print("\n[Step 2] Loading terminal sequences...")

    try:
        uvig_ends = pd.read_csv(
            config['data']['terminal_sequences'],
            sep='\t',
            header=0,
            names=['full_name', 'start_25bp', 'end_25bp']
        )
        uvig_ends['uvig'] = uvig_ends['full_name'].str.split('|').str[0]
        print(f"✓ Loaded terminal sequences for {len(uvig_ends):,} UViGs")
    except FileNotFoundError:
        print("⚠ Terminal sequence file not found - using mock data for demo")
        # Create mock data for demonstration
        uvig_ends = pd.DataFrame({
            'full_name': df_uvigs['uvig'].head(100),
            'start_25bp': ['ACGTACGTACGTACGTACGTACGTA'] * 100,
            'end_25bp': ['TACGTACGTACGTACGTACGTACGT'] * 100,
        })
        uvig_ends['uvig'] = uvig_ends['full_name']
        print(f"  Using mock data: {len(uvig_ends)} sequences")

    # ====================================================================
    # Filter terminal sequences to our UViGs
    # ====================================================================
    print("\n[Step 3] Filtering terminal sequences to our dataset...")

    uvig_ends_filtered = uvig_ends[uvig_ends['uvig'].isin(df_uvigs['uvig'].unique())]
    print(f"✓ Found terminal sequences for {len(uvig_ends_filtered):,} UViGs")

    # ====================================================================
    # Flag representative UViGs
    # ====================================================================
    print("\n[Step 4] Flagging representative UViGs...")

    df_uvigs = flag_representative_uvigs(df_uvigs, representative)
    n_rep = df_uvigs['is_representative'].sum()
    print(f"✓ Identified {n_rep:,} representative UViGs")

    # ====================================================================
    # Compute terminal sequence mismatches
    # ====================================================================
    print("\n[Step 5] Computing terminal sequence mismatches...")

    mismatch_results = compute_uvig_mismatches_per_member(
        uvig_ends_filtered, df_uvigs,
        start_col=config['terminal_analysis']['start_column'],
        end_col=config['terminal_analysis']['end_column'],
        max_mismatches=config['terminal_analysis']['max_mismatches']
    )

    print(f"✓ Analyzed {len(mismatch_results):,} UViGs")

    # Show distribution of mismatches
    print("\n  Mismatch distribution (start positions):")
    print(mismatch_results['mismatches_start'].describe())

    print("\n  Mismatch distribution (end positions):")
    print(mismatch_results['mismatches_end'].describe())

    # ====================================================================
    # Filter by passing members
    # ====================================================================
    print("\n[Step 6] Filtering vOTUs by passing member count...")

    min_samples = config['filtering']['min_unique_samples']

    summary, filtered_v2 = filter_by_passing_members(
        mismatch_results,
        min_unique_samples=min_samples
    )

    print(f"✓ Filtered to vOTUs with ≥{min_samples} passing members")
    print(f"  Result: {filtered_v2['votu'].nunique():,} passing vOTUs")

    # ====================================================================
    # Summarize by thresholds
    # ====================================================================
    print("\n[Step 7] Summary statistics by sample threshold...")

    summary_stats = summarize_passing_uvigs(summary, thresholds=[2, 3, 4, 5])
    print("\n  Counts at different thresholds:")
    for col in summary_stats.columns:
        count = summary_stats[col].values[0]
        print(f"    {col:15s}: {count:5,} vOTUs")

    # ====================================================================
    # Quality assessment
    # ====================================================================
    print("\n[Step 8] Quality assessment of passing sequences...")

    passing_by_quality = filtered_v2.groupby('quality')['votu'].nunique()
    print("\n  Passing vOTUs by quality:")
    for quality, count in passing_by_quality.items():
        print(f"    {quality:20s}: {count:5,}")

    passing_by_source = filtered_v2.groupby('source')['votu'].nunique()
    print("\n  Passing vOTUs by source:")
    for source, count in passing_by_source.items():
        print(f"    {source:20s}: {count:5,}")

    # ====================================================================
    # Conservation patterns
    # ====================================================================
    print("\n[Step 9] Terminal sequence conservation patterns...")

    both_pass = (
        mismatch_results['both_pass'] &
        (~mismatch_results['is_representative'])
    ).sum()

    start_only = (
        ((mismatch_results['mismatches_start'] <= config['terminal_analysis']['max_mismatches']) &
         (mismatch_results['mismatches_end'] > config['terminal_analysis']['max_mismatches'])) &
        (~mismatch_results['is_representative'])
    ).sum()

    end_only = (
        ((mismatch_results['mismatches_start'] > config['terminal_analysis']['max_mismatches']) &
         (mismatch_results['mismatches_end'] <= config['terminal_analysis']['max_mismatches'])) &
        (~mismatch_results['is_representative'])
    ).sum()

    total_nonrep = (~mismatch_results['is_representative']).sum()

    print(f"\n  Among {total_nonrep:,} non-representative UViGs:")
    print(f"    Both termini conserved: {both_pass:7,} ({both_pass/total_nonrep*100:5.1f}%)")
    print(f"    Start only conserved:   {start_only:7,} ({start_only/total_nonrep*100:5.1f}%)")
    print(f"    End only conserved:     {end_only:7,} ({end_only/total_nonrep*100:5.1f}%)")

    # ====================================================================
    # Save results
    # ====================================================================
    print("\n[Step 10] Saving results...")

    output_dir = PROJECT_ROOT / 'results'
    output_dir.mkdir(exist_ok=True)

    # Save mismatch analysis
    mismatch_output = output_dir / 'terminal_mismatch_analysis.tsv'
    mismatch_results.to_csv(mismatch_output, sep='\t', index=False)
    print(f"✓ Saved mismatch analysis to {mismatch_output}")

    # Save passing vOTUs
    passing_output = output_dir / 'passing_terminal_sequences.tsv'
    filtered_v2.to_csv(passing_output, sep='\t', index=False)
    print(f"✓ Saved passing sequences to {passing_output}")

    # Save summary
    summary_output = output_dir / 'terminal_summary.tsv'
    summary.to_csv(summary_output, sep='\t', index=False)
    print(f"✓ Saved summary to {summary_output}")

    # ====================================================================
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
