"""
Filtering module for viral sequence datasets.

This module contains vectorized filtering functions for:
- Length-based filtering (±bp tolerance)
- Group size filtering (minimum UViGs per vOTU)
- Quality/taxonomy-based filtering
"""

import pandas as pd
import numpy as np
from typing import Dict, Tuple, List


def filter_exclusive_linear_votus(df: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
    """
    Identify vOTUs that contain ONLY linear topology UViGs.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata DataFrame with 'votu' and 'uvig_topology' columns.

    Returns
    -------
    tuple
        (filtered_df, n_exclusive_votu)
        - filtered_df: DataFrame with linear-only vOTUs
        - n_exclusive_votu: Count of exclusively linear vOTUs
    """
    # Count unique topologies per vOTU
    exclusively_linear_votus = (
        df.groupby("votu")["uvig_topology"]
        .nunique()
        .loc[lambda x: (x == 1)]
        .index
    )

    # Filter to those with Linear topology
    linear_votus = df[df["votu"].isin(exclusively_linear_votus)]
    linear_votus = linear_votus[linear_votus["uvig_topology"] == "Linear"]

    exclusive_linear_votus = linear_votus["votu"].unique()

    return linear_votus, len(exclusive_linear_votus)


def remove_nested_votus(
    df: pd.DataFrame,
    exclude_uvig_list: pd.DataFrame
) -> Tuple[pd.DataFrame, int]:
    """
    Remove vOTUs that contain fragmented/shorter UViGs falling into longer vOTUs.

    Parameters
    ----------
    df : pd.DataFrame
        Linear vOTU DataFrame.
    exclude_uvig_list : pd.DataFrame
        DataFrame with 'uvig' column containing UViGs to exclude.

    Returns
    -------
    tuple
        (cleaned_df, n_removed_votus)
    """
    nested_votus = df[
        df['uvig'].isin(exclude_uvig_list['uvig'].to_list())
    ]['votu'].unique()

    cleaned = df[~df['votu'].isin(nested_votus)].copy()
    n_removed = len(nested_votus)

    return cleaned, n_removed


def filter_by_length_difference(
    df: pd.DataFrame,
    correspondence: pd.DataFrame,
    max_bp_diff: int = 25
) -> pd.DataFrame:
    """
    Filter UViGs based on length difference from representative.

    Keep only UViGs whose length is within ±max_bp_diff of their representative.

    Parameters
    ----------
    df : pd.DataFrame
        Linear UViG DataFrame with 'uvig', 'votu', 'length' columns.
    correspondence : pd.DataFrame
        Mapping of vOTU to representative uvig with 'votu', 'uvig', 'length' columns.
    max_bp_diff : int, optional
        Maximum allowed absolute length difference in bp (default: 25).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only length-similar UViGs.
    """
    # Create representative length mapping
    corr_subset = correspondence[['votu', 'uvig', 'length']].rename(
        columns={'uvig': 'rep_uvig', 'length': 'rep_length'}
    )

    # Merge representative info
    df_merged = df.merge(corr_subset, on='votu', how='left')

    # Vectorized length filtering
    mask_length = (df_merged["length"] - df_merged["rep_length"]).abs() <= max_bp_diff

    return df_merged[mask_length][df.columns].copy()


def filter_by_group_size(
    df: pd.DataFrame,
    min_uvigs_per_votu: int = 3
) -> pd.DataFrame:
    """
    Remove vOTUs that have fewer than min_uvigs_per_votu members.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with 'votu' column.
    min_uvigs_per_votu : int, optional
        Minimum number of UViGs per vOTU (default: 3).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only sufficiently large vOTU groups.
    """
    votu_counts = df.groupby("votu")["uvig"].transform("count")
    return df[votu_counts >= min_uvigs_per_votu].copy().reset_index(drop=True)


def filter_linear_uvigs_comprehensive(
    linear_votus: pd.DataFrame,
    quality_tier_df: pd.DataFrame,
    correspondence: pd.DataFrame,
    max_bp_diff: int = 25,
    min_uvigs_per_votu: int = 3
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """
    Comprehensive filtering pipeline for linear UViGs.

    This is the main filtering function combining all individual filters.

    Parameters
    ----------
    linear_votus : pd.DataFrame
        All exclusively linear vOTUs.
    quality_tier_df : pd.DataFrame
        Quality tier subset (e.g., High-quality only).
    correspondence : pd.DataFrame
        vOTU to representative UViG mapping.
    max_bp_diff : int, optional
        Max bp difference from representative (default: 25).
    min_uvigs_per_votu : int, optional
        Min UViGs per vOTU (default: 3).

    Returns
    -------
    tuple
        (filtered_df, stats_dict)
        stats_dict contains counts at each filtering step
    """
    stats = {}

    # Step 1: Filter to quality tier only
    df = linear_votus[linear_votus["votu"].isin(quality_tier_df["votu"])].copy()
    stats['after_quality_filter'] = {
        'votu': df['votu'].nunique(),
        'uvig': df['uvig'].nunique()
    }

    # Step 2: Filter by length
    df = filter_by_length_difference(df, correspondence, max_bp_diff)
    stats['after_length_filter'] = {
        'votu': df['votu'].nunique(),
        'uvig': df['uvig'].nunique()
    }

    # Step 3: Filter by group size
    df = filter_by_group_size(df, min_uvigs_per_votu)
    stats['after_group_size_filter'] = {
        'votu': df['votu'].nunique(),
        'uvig': df['uvig'].nunique()
    }

    return df, stats


def filter_by_min_genome_length(
    df: pd.DataFrame,
    min_length_dict: Dict[str, int],
    class_col: str = "class",
    length_col: str = "length"
) -> pd.DataFrame:
    """
    Filter genomes by class-specific minimum length.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with taxonomy class and genome length columns.
    min_length_dict : dict
        Mapping of class name -> minimum length in bp.
        Example: {'Caudoviricetes': 12275, 'Leviviricetes': 3000}
    class_col : str, optional
        Name of class column (default: 'class').
    length_col : str, optional
        Name of length column (default: 'length').

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only genomes meeting minimum length requirements.
    """
    min_length_series = df[class_col].map(min_length_dict).fillna(0)
    return df.loc[df[length_col] >= min_length_series].copy()


def identify_length_outliers(
    df: pd.DataFrame,
    votu_col: str = "votu",
    length_col: str = "length",
    max_delta_bp: int = 1000
) -> pd.DataFrame:
    """
    Identify vOTUs with large length variations (potential assemblies from different species).

    Parameters
    ----------
    df : pd.DataFrame
        Full metadata file with all UViGs.
    votu_col : str, optional
        vOTU column name (default: 'votu').
    length_col : str, optional
        Length column name (default: 'length').
    max_delta_bp : int, optional
        Maximum allowed difference between longest and shortest UViG (default: 1000).

    Returns
    -------
    pd.DataFrame
        Summary of vOTUs with length > max_delta_bp, sorted by largest difference.
    """
    summary = (
        df.sort_values(by=[votu_col, length_col])
        .groupby(votu_col)
        .agg({
            'uvig': ['first', 'last'],
            length_col: ['min', 'max']
        })
        .reset_index()
    )

    summary.columns = ['votu', 'shortest_uvig', 'longest_uvig', 'min_length', 'max_length']
    summary['length_delta'] = summary['max_length'] - summary['min_length']

    return summary[summary['length_delta'] > max_delta_bp].sort_values(
        by='length_delta', ascending=False
    )


def print_filter_summary(stats: Dict) -> None:
    """
    Print a formatted summary of filtering statistics.

    Parameters
    ----------
    stats : dict
        Statistics dictionary from filter_linear_uvigs_comprehensive.
    """
    print("\n" + "="*60)
    print("FILTERING SUMMARY")
    print("="*60)

    for step, counts in stats.items():
        print(f"\n{step}:")
        print(f"  vOTUs: {counts['votu']:,}")
        print(f"  UViGs: {counts['uvig']:,}")
