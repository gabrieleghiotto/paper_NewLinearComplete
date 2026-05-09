"""
Terminal sequence analysis module.

Analyzes conserved start/end sequences (termini) of viral genomes to identify
orthologous sequences across members of the same vOTU.
"""

import pandas as pd
import numpy as np
from typing import Tuple, Dict


def flag_representative_uvigs(
    final_set: pd.DataFrame,
    correspondence: pd.DataFrame,
    uvig_col: str = 'uvig',
    votu_col: str = 'votu'
) -> pd.DataFrame:
    """
    Mark which UViGs are representatives of their vOTU.

    Parameters
    ----------
    final_set : pd.DataFrame
        UViG dataset with 'uvig' and 'votu' columns.
    correspondence : pd.DataFrame
        Representative mapping with 'votu', 'uvig', 'length' columns.
    uvig_col : str, optional
        Name of UViG ID column (default: 'uvig').
    votu_col : str, optional
        Name of vOTU ID column (default: 'votu').

    Returns
    -------
    pd.DataFrame
        Input dataframe with added 'is_representative' boolean column.
    """
    df = final_set.copy()

    # Create mapping
    rep_dict = dict(zip(correspondence[votu_col], correspondence[uvig_col]))

    # Vectorized comparison
    df['is_representative'] = df.apply(
        lambda row: row[uvig_col] == rep_dict.get(row[votu_col], None),
        axis=1
    )

    return df


def compute_sequence_mismatches(
    query_seq: str,
    ref_seq: str,
    min_length: int | None = None
) -> int:
    """
    Count mismatches between two sequences of equal length.

    Parameters
    ----------
    query_seq : str
        Query DNA sequence.
    ref_seq : str
        Reference DNA sequence.
    min_length : int, optional
        If sequences differ in length, trim to this length before comparison.

    Returns
    -------
    int
        Number of mismatches.
    """
    if pd.isna(query_seq) or pd.isna(ref_seq):
        return -1

    query_seq = str(query_seq).strip()
    ref_seq = str(ref_seq).strip()

    if not query_seq or not ref_seq:
        return -1

    # Trim to equal length
    min_len = min(len(query_seq), len(ref_seq))
    if min_length:
        min_len = min(min_len, min_length)

    query_array = np.array(list(query_seq[:min_len]))
    ref_array = np.array(list(ref_seq[:min_len]))

    return (query_array != ref_array).sum()


def compute_uvig_mismatches_per_member(
    uvig_ends: pd.DataFrame,
    final_set: pd.DataFrame,
    start_col: str,
    end_col: str,
    max_mismatches: int = 1
) -> pd.DataFrame:
    """
    Compute terminal sequence mismatches for each UViG against its representative.

    Parameters
    ----------
    uvig_ends : pd.DataFrame
        DataFrame with 'uvig' and sequence columns (start_col, end_col).
    final_set : pd.DataFrame
        UViGs to analyze with 'uvig', 'votu', 'is_representative' columns.
    start_col : str
        Name of start sequence column in uvig_ends.
    end_col : str
        Name of end sequence column in uvig_ends.
    max_mismatches : int, optional
        Threshold for determining pass/fail (default: 1).

    Returns
    -------
    pd.DataFrame
        Analysis table with columns:
        ['uvig', 'votu', 'seq_start', 'seq_end', 'mismatches_start',
         'mismatches_end', 'is_representative', 'both_pass', 'taxon_oid']
    """
    # Prepare data
    df_set = final_set.copy()
    df_ends = uvig_ends.copy()

    # Clean UViG names
    def clean_uvig(name):
        if isinstance(name, str):
            return name.strip().lstrip('>').split('|')[0].replace('.fasta', '')
        return str(name)

    df_set['uvig'] = df_set['uvig'].apply(clean_uvig)
    df_ends['uvig'] = df_ends['uvig'].apply(clean_uvig)

    # Get representatives and their sequences
    reps = df_set[df_set['is_representative']][['votu', 'uvig']].copy()
    reps = reps.merge(df_ends[['uvig', start_col, end_col]], on='uvig', how='left')
    reps = reps.rename(columns={start_col: 'rep_start', end_col: 'rep_end'})

    # Merge representative sequences into all members
    merged = df_set.merge(reps[['votu', 'rep_start', 'rep_end']], on='votu', how='left')
    merged = merged.merge(df_ends[['uvig', start_col, end_col]], on='uvig', how='left')
    merged = merged.rename(columns={start_col: 'seq_start', end_col: 'seq_end'})

    # Fill NaN with empty strings
    seq_cols = ['seq_start', 'seq_end', 'rep_start', 'rep_end']
    for col in seq_cols:
        merged[col] = merged[col].fillna('').astype(str)

    # Process non-representatives
    mask_nonrep = ~merged['is_representative']
    merged_nonrep = merged[mask_nonrep].copy()

    # Skip invalid rows
    valid_mask = (
        (merged_nonrep['seq_start'] != '') &
        (merged_nonrep['seq_end'] != '') &
        (merged_nonrep['rep_start'] != '') &
        (merged_nonrep['rep_end'] != '')
    )
    merged_nonrep = merged_nonrep[valid_mask].copy()

    # Compute mismatches vectorized
    start_array = np.array([list(s) for s in merged_nonrep['seq_start']])
    end_array = np.array([list(s) for s in merged_nonrep['seq_end']])
    rep_start_array = np.array([list(s) for s in merged_nonrep['rep_start']])
    rep_end_array = np.array([list(s) for s in merged_nonrep['rep_end']])

    # Align lengths
    min_start = min(start_array.shape[1], rep_start_array.shape[1])
    min_end = min(end_array.shape[1], rep_end_array.shape[1])
    start_array = start_array[:, :min_start]
    rep_start_array = rep_start_array[:, :min_start]
    end_array = end_array[:, :min_end]
    rep_end_array = rep_end_array[:, :min_end]

    merged_nonrep['mismatches_start'] = (start_array != rep_start_array).sum(axis=1)
    merged_nonrep['mismatches_end'] = (end_array != rep_end_array).sum(axis=1)
    merged_nonrep['both_pass'] = (
        (merged_nonrep['mismatches_start'] <= max_mismatches) &
        (merged_nonrep['mismatches_end'] <= max_mismatches)
    )

    # Representatives automatically pass
    reps_only = merged[merged['is_representative']].copy()
    reps_only['mismatches_start'] = 0
    reps_only['mismatches_end'] = 0
    reps_only['both_pass'] = True

    # Combine
    final_table = pd.concat([merged_nonrep, reps_only], ignore_index=True)
    final_table = final_table[[
        'uvig', 'votu', 'seq_start', 'seq_end',
        'mismatches_start', 'mismatches_end',
        'is_representative', 'both_pass', 'taxon_oid'
    ]]

    return final_table


def filter_by_passing_members(
    final_table: pd.DataFrame,
    min_unique_samples: int = 3
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter vOTUs based on number of passing non-representative members from different samples.

    Parameters
    ----------
    final_table : pd.DataFrame
        Mismatch analysis table from compute_uvig_mismatches_per_member.
    min_unique_samples : int, optional
        Minimum number of unique samples with passing sequences (default: 3).

    Returns
    -------
    tuple
        (summary_df, filtered_table)
        - summary_df: Count of passing unique samples per vOTU
        - filtered_table: Original table filtered to qualifying vOTUs
    """
    # Count passing members per vOTU from unique samples
    mask = (~final_table['is_representative']) & (final_table['both_pass'] == True)
    passing_nonrep = final_table[mask].copy()

    summary_per_votu = (
        passing_nonrep.groupby('votu')['taxon_oid']
        .nunique()
        .reset_index()
        .rename(columns={'taxon_oid': 'n_unique_samples_passing'})
    )

    # Include all vOTUs (fill 0 for those with no passing members)
    all_votus = pd.DataFrame({'votu': final_table['votu'].unique()})
    summary_per_votu = all_votus.merge(summary_per_votu, on='votu', how='left').fillna(0)
    summary_per_votu['n_unique_samples_passing'] = summary_per_votu[
        'n_unique_samples_passing'
    ].astype(int)

    # Filter
    votus_validi = summary_per_votu[
        summary_per_votu['n_unique_samples_passing'] >= min_unique_samples
    ]['votu']

    filtered_table = final_table[final_table['votu'].isin(votus_validi)].copy()

    return summary_per_votu, filtered_table


def summarize_passing_uvigs(summary_df: pd.DataFrame, thresholds: list = None) -> pd.DataFrame:
    """
    Generate summary statistics for passing UViGs at different thresholds.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Summary output from filter_by_passing_members.
    thresholds : list, optional
        Sample count thresholds to report (default: [2, 3, 4, 5]).

    Returns
    -------
    pd.DataFrame
        Summary with counts of vOTUs meeting each threshold.
    """
    if thresholds is None:
        thresholds = [2, 3, 4, 5]

    summary_dict = {
        f'>={t}_samples': (summary_df['n_unique_samples_passing'] >= t).sum()
        for t in thresholds
    }

    return pd.DataFrame([summary_dict])
