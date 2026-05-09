"""
Taxonomy extraction and classification module.

Handles parsing ICTV and host taxonomy strings, extracting taxonomic ranks,
and classifying genomes by type.
"""

import pandas as pd
from typing import Dict, Tuple


def extract_taxonomy_level(tax_string: str, prefix: str) -> str:
    """
    Extract a specific taxonomic level from a semicolon-delimited taxonomy string.

    Parameters
    ----------
    tax_string : str
        Full taxonomy string (e.g., "d__Bacteria;p__Bacillota;c__...")
    prefix : str
        Taxonomic prefix to find (e.g., 'g__' for genus, 's__' for species).

    Returns
    -------
    str
        Extracted taxonomy level with prefix removed, or "Unknown" if not found.

    Examples
    --------
    >>> tax = "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus"
    >>> extract_taxonomy_level(tax, 'g__')
    'Bacillus'
    """
    if pd.isna(tax_string) or tax_string == "NaN":
        return "Unknown"

    ranks = str(tax_string).split(';')
    for rank in ranks:
        if rank.startswith(prefix):
            return rank.replace(prefix, '').replace('_', ' ').strip()

    return "Not assigned"


def extract_ictv_and_host_class(
    df: pd.DataFrame,
    ictv_col: str = "ictv_taxonomy",
    host_col: str = "host_taxonomy",
    class_fill: str = "c__",
    host_class_fill: str = "Unclassified"
) -> pd.DataFrame:
    """
    Extract viral and host class from ICTV and host taxonomy strings.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with taxonomy columns.
    ictv_col : str, optional
        Column containing ICTV taxonomy (default: 'ictv_taxonomy').
    host_col : str, optional
        Column containing host taxonomy (default: 'host_taxonomy').
    class_fill : str, optional
        Value for missing viral class (default: 'c__').
    host_class_fill : str, optional
        Value for missing host class (default: 'Unclassified').

    Returns
    -------
    pd.DataFrame
        Copy of df with added 'class' and 'host_class' columns.
    """
    df = df.copy()

    # Extract viral class
    df["class"] = (
        df[ictv_col]
        .str.extract(r"c__([^;]*)", expand=False)
        .replace("", pd.NA)
        .fillna(class_fill)
    )

    # Extract host class
    df["host_class"] = (
        df[host_col]
        .str.extract(r"c__([^;]*)", expand=False)
        .replace("", pd.NA)
        .fillna(host_class_fill)
    )

    # Mark unclassified viral genomes
    df.loc[df[ictv_col].isna(), "class"] = "unclassified"

    return df


def extract_order_family(df: pd.DataFrame, col: str = "ictv_taxonomy") -> pd.DataFrame:
    """
    Extract taxonomic order, family, and genus from ICTV taxonomy string.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with ICTV taxonomy column.
    col : str, optional
        Name of taxonomy column (default: 'ictv_taxonomy').

    Returns
    -------
    pd.DataFrame
        Copy of df with added 'order', 'family', 'genus' columns.

    Examples
    --------
    >>> df = extract_order_family(df)
    >>> df[['order', 'family', 'genus']]
    """
    def parse_tax(tax):
        order, family, genus = None, None, None

        if pd.isna(tax):
            return pd.Series([order, family, genus])

        for level in str(tax).split(";"):
            if level.startswith("o__"):
                order = level.replace("o__", "").replace("_", " ").strip()
            elif level.startswith("f__"):
                family = level.replace("f__", "").replace("_", " ").strip()
            elif level.startswith("g__"):
                genus = level.replace("g__", "").replace("_", " ").strip()

        return pd.Series([order, family, genus])

    df = df.copy()
    df[["order", "family", "genus"]] = df[col].apply(parse_tax)

    return df


def classify_genome_type(viral_class: str) -> str:
    """
    Classify viral genome as DNA, RNA, or Unknown based on ICTV class.

    Parameters
    ----------
    viral_class : str
        ICTV viral class (classiviricetes name).

    Returns
    -------
    str
        'DNA', 'RNA', or 'Unknown'.

    References
    ----------
    Based on ICTV classification from the Virus Metadata Resource (ICTV 2023).
    """
    rna_classes = {
        "Ainoaviricetes", "Alsuviricetes", "Amabiliviricetes",
        "Bunyaviricetes", "Chrymotiviricetes", "Duplopiviricetes",
        "Howeltoviricetes", "Leviviricetes", "Magsaviricetes",
        "Miaviricetes", "Milneviricetes", "Monjiviricetes",
        "Pisoniviricetes", "Resentoviricetes", "Revtraviricetes",
        "Stelpaviricetes", "Vidaverviricetes", "Tolucaviricetes",
    }

    dna_classes = {
        "Arfiviricetes", "Caudoviricetes", "Faserviricetes",
        "Herviviricetes", "Huolimaviricetes", "Laserviricetes",
        "Malgrandaviricetes", "Maveriviricetes", "Megaviricetes",
        "Papovaviricetes", "Polintoviricetes", "Pokkesviricetes",
        "Quintoviricetes", "Repensiviricetes", "Tectiliviricetes",
        "Tokiviricetes"
    }

    if viral_class in rna_classes:
        return "RNA"
    elif viral_class in dna_classes:
        return "DNA"
    else:
        return "Unknown"


def add_genome_type_column(df: pd.DataFrame, class_col: str = "class") -> pd.DataFrame:
    """
    Add genome_type column to dataframe based on viral class.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with viral class column.
    class_col : str, optional
        Name of class column (default: 'class').

    Returns
    -------
    pd.DataFrame
        Copy of df with added 'Genome_type' column.
    """
    df = df.copy()
    df['Genome_type'] = df[class_col].apply(classify_genome_type)
    return df


def get_min_genome_lengths() -> Dict[str, int]:
    """
    Return dictionary of minimum genome lengths by ICTV class.

    Based on class-specific length requirements from viral databases.

    Returns
    -------
    dict
        Mapping of class name -> minimum genome length in bp.
    """
    return {
        "Alsuviricetes": 1000,
        "Amabiliviricetes": 2000,
        "Arfiviricetes": 1000,
        "Bunyaviricetes": 1700,
        "Caudoviricetes": 12275,
        "Chrymotiviricetes": 2700,
        "Duplopiviricetes": 1400,
        "Faserviricetes": 4500,
        "Howeltoviricetes": 2400,
        "Huolimaviricetes": 1000,
        "Laserviricetes": 16000,
        "Leviviricetes": 3000,
        "Magsaviricetes": 400,
        "Malgrandaviricetes": 4400,
        "Maveriviricetes": 14500,
        "Megaviricetes": 70000,
        "Miaviricetes": 900,
        "Milneviricetes": 1400,
        "Monjiviricetes": 8500,
        "Naldaviricetes": 80000,
        "Papovaviricetes": 5000,
        "Pisoniviricetes": 900,
        "Pokkesviricetes": 128000,
        "Polintoviricetes": 1200,
        "Quintoviricetes": 4000,
        "Repensiviricetes": 2000,
        "Resentoviricetes": 200,
        "Revtraviricetes": 3200,
        "Stelpaviricetes": 2300,
        "Tectiliviricetes": 8700,
        "Tokiviricetes": 15900,
        "Tolucaviricetes": 1400,
        "Vidaverviricetes": 2300,
    }


def get_class_summary(
    df: pd.DataFrame,
    class_col: str = "class"
) -> pd.Series:
    """
    Get count of vOTUs per viral class.

    Parameters
    ----------
    df : pd.DataFrame
        Master dataframe with vOTU and class information.
    class_col : str, optional
        Name of class column (default: 'class').

    Returns
    -------
    pd.Series
        Counts of unique vOTUs per class, sorted descending.
    """
    return df.groupby(class_col)['votu'].nunique().sort_values(ascending=False)
