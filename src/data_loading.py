"""
Data loading and basic metadata extraction module.

This module handles loading viral metadata from TSV files and performing
initial filtering and data organization.
"""

import pandas as pd
from pathlib import Path
from typing import Tuple, Dict, Any


def load_metadata(metadata_path: str | Path, sep: str = '\t') -> pd.DataFrame:
    """
    Load metaVR metadata file.

    Parameters
    ----------
    metadata_path : str or Path
        Path to UVIG metadata TSV file.
    sep : str, optional
        Separator used in the file (default: '\t').

    Returns
    -------
    pd.DataFrame
        Loaded metadata with columns: uvig, votu, length, uvig_topology, etc.
    """
    df = pd.read_csv(metadata_path, sep=sep)
    return df


def load_representative_uvigs(header_path: str | Path) -> pd.DataFrame:
    """
    Load representative UViG identifiers.

    Representative UViGs are the canonical sequence for each vOTU.

    Parameters
    ----------
    header_path : str or Path
        Path to header file containing representative UViG names.

    Returns
    -------
    pd.DataFrame
        DataFrame with single column 'uvig' containing representative UViG names.
    """
    df = pd.read_csv(header_path, header=None, names=['uvig'], dtype=str)
    return df


def load_excluded_sequences(exclude_path: str | Path) -> pd.DataFrame:
    """
    Load sequences that should be excluded from analysis (e.g., fragmented sequences).

    Parameters
    ----------
    exclude_path : str or Path
        Path to file containing UViG names to exclude.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns 'name' and 'uvig'.
    """
    df = pd.read_csv(exclude_path, header=None, names=['name'], dtype=str)
    df['uvig'] = df['name'].str.split('|').str[0]
    return df


def extract_topology_quality_subsets(
    df: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """
    Split metadata into subsets based on topology and quality.

    Parameters
    ----------
    df : pd.DataFrame
        Full metadata DataFrame with 'uvig_topology' and 'quality' columns.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'circular_complete': direct or inverted terminal repeats, Complete quality
        - 'linear_complete': Linear topology, Complete quality
        - 'all_linear': All linear topology UViGs
    """
    subsets = {}

    # Circular complete (any terminal repeat type)
    circular_mask = (
        (df['uvig_topology'] == 'Direct terminal repeat') |
        (df['uvig_topology'] == 'Inverted terminal repeat')
    )
    subsets['circular_complete'] = df[circular_mask & (df['quality'] == 'Complete')].copy()

    # Linear complete
    linear_mask = df['uvig_topology'] == 'Linear'
    subsets['linear_complete'] = df[linear_mask & (df['quality'] == 'Complete')].copy()

    # All linear
    subsets['all_linear'] = df[linear_mask].copy()

    return subsets


def extract_quality_tiers(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Extract quality tier subsets from linear genomes.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata DataFrame with 'quality' and 'viral_confidence' columns.
        Typically the result of filtering all_linear UViGs.

    Returns
    -------
    dict
        Dictionary with keys:
        - 'high_quality': quality == 'High-quality'
        - 'medium_quality': quality == 'Medium-quality'
        - 'not_determined_high': quality == 'Not-determined' AND viral_confidence == 'High'
        - 'not_determined_low': quality == 'Not-determined' AND viral_confidence == 'Low'
    """
    subsets = {}

    subsets['high_quality'] = df[df['quality'] == 'High-quality'].copy()
    subsets['medium_quality'] = df[df['quality'] == 'Medium-quality'].copy()

    nd_mask = df['quality'] == 'Not-determined'
    subsets['not_determined_high'] = df[nd_mask & (df['viral_confidence'] == 'High')].copy()
    subsets['not_determined_low'] = df[nd_mask & (df['viral_confidence'] == 'Low')].copy()

    return subsets


def get_votu_representatives(
    df: pd.DataFrame,
    representative_uvigs: pd.DataFrame
) -> pd.DataFrame:
    """
    Filter to only representative UViGs from the full dataset.

    Parameters
    ----------
    df : pd.DataFrame
        Full metadata DataFrame.
    representative_uvigs : pd.DataFrame
        DataFrame with 'uvig' column containing representative UViG names.

    Returns
    -------
    pd.DataFrame
        Subset of df containing only representative UViGs.
    """
    return df[df['uvig'].isin(representative_uvigs['uvig'].to_list())].copy()


def count_votus_per_category(
    df: pd.DataFrame,
    group_col: str = 'source'
) -> pd.Series:
    """
    Count unique vOTUs per category.

    Parameters
    ----------
    df : pd.DataFrame
        Metadata DataFrame with 'votu' column.
    group_col : str, optional
        Column to group by (default: 'source').

    Returns
    -------
    pd.Series
        Count of unique vOTUs per category.
    """
    return df.groupby(group_col)['votu'].nunique()


def save_results(
    df: pd.DataFrame,
    output_path: str | Path,
    sep: str = '\t'
) -> None:
    """
    Save filtered results to disk.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to save.
    output_path : str or Path
        Output file path.
    sep : str, optional
        Separator to use (default: '\t').
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep=sep, index=False)
    print(f"✓ Saved {len(df)} rows to {output_path}")
