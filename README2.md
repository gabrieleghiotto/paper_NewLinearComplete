# Viral Linear Complete Genomes Analysis

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A reproducible pipeline for analyzing viral genomes with linear topology, including filtering, quality assessment, taxonomy classification, and terminal sequence conservation analysis.

## Overview

This project provides tools to:

- **Filter** viral sequences by topology, quality, and genome length
- **Identify** conserved terminal sequences across genome variants
- **Classify** genomes by taxonomy (ICTV classes) and genome type (DNA/RNA)
- **Generate** publication-quality figures and summary statistics

The pipeline is designed for the metaVR dataset but is adaptable to any viral metadata in TSV format.

## Features

✨ **Modular Architecture** - Functions split into logical modules for reusability
✨ **Configuration-Driven** - YAML config files instead of hardcoded parameters
✨ **Reproducible** - Deterministic results with documented algorithms
✨ **Well-Documented** - Docstrings and type hints throughout
✨ **Vectorized** - Fast pandas/numpy operations for large datasets

## Installation

### Requirements
- Python 3.8+
- ~2 GB RAM (for full metaVR dataset)

### Setup

```bash
# Clone repository
git clone https://github.com/yourusername/viral-complete-genomes.git
cd viral-complete-genomes

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### 1. Configure your environment

Edit `config/config.yaml` with your data paths:

```yaml
data:
  input_metadata: "/path/to/UVIG_metadata.tsv"
  header_repr: "/path/to/representatives.txt"
  # ... other paths
```

### 2. Run the pipeline

```bash
# Run complete analysis
python scripts/run_pipeline.py --config config/config.yaml

# Or specify custom output directory
python scripts/run_pipeline.py --config config/config.yaml --output-dir ./my_results
```

### 3. View results

```
results/
├── tables/
│   ├── final_filtered_linear_complete.tsv
│   ├── linear_complete_high_quality.tsv
│   └── summary_statistics.tsv
└── figures/
    ├── Fig2_overview.svg
    ├── Fig4_5_taxonomy.svg
    └── ...
```

## Data Requirements

Input TSV file should contain at minimum:

| Column | Description |
|--------|-------------|
| `uvig` | Unique UViG identifier |
| `votu` | Viral operational taxonomic unit |
| `length` | Genome length in bp |
| `uvig_topology` | Genome topology (Linear, Direct terminal repeat, Inverted terminal repeat) |
| `quality` | Quality assessment (High-quality, Medium-quality, Not-determined) |
| `ictv_taxonomy` | ICTV taxonomy string (semicolon-delimited) |
| `host_taxonomy` | Host taxonomy string (semicolon-delimited) |
| `source` | Data source (Metagenome, Metatranscriptome, Isolate, RefSeq) |

## Module Documentation

### `src/data_loading.py`
Load and organize viral metadata from multiple sources.

```python
from src.data_loading import load_metadata, extract_topology_quality_subsets

df = load_metadata('data.tsv')
subsets = extract_topology_quality_subsets(df)
circular_complete = subsets['circular_complete']
```

### `src/filtering.py`
Comprehensive filtering pipeline with multiple strategies.

```python
from src.filtering import filter_linear_uvigs_comprehensive

filtered, stats = filter_linear_uvigs_comprehensive(
    linear_votus, quality_tier, correspondence,
    max_bp_diff=25, min_uvigs_per_votu=3
)
print_filter_summary(stats)
```

### `src/taxonomy.py`
Extract and classify viral taxonomy.

```python
from src.taxonomy import extract_ictv_and_host_class, add_genome_type_column

df = extract_ictv_and_host_class(df)  # Extracts class from taxonomy
df = add_genome_type_column(df)       # Classifies as DNA/RNA/Unknown
```

### `src/terminal_analysis.py`
Analyze conserved terminal sequences.

```python
from src.terminal_analysis import compute_uvig_mismatches_per_member

mismatches = compute_uvig_mismatches_per_member(
    uvig_ends, final_set,
    start_col='start_25bp', end_col='end_25bp',
    max_mismatches=3
)
```

## Configuration

Key configuration sections in `config/config.yaml`:

**filtering**
```yaml
filtering:
  max_bp_difference: 25          # ±bp from representative
  min_uvigs_per_votu: 3          # Minimum group size
  min_unique_samples: 3          # For terminal analysis
```

**taxonomy**
```yaml
taxonomy:
  min_length_dict:
    Caudoviricetes: 12275
    Leviviricetes: 3000
    # ... etc
```

See `config/config.yaml` for full configuration documentation.

## Example: Custom Analysis

```python
import pandas as pd
from pathlib import Path
import yaml

# Load configuration
with open('config/config.yaml') as f:
    config = yaml.safe_load(f)

# Import modules
from src.data_loading import load_metadata, get_votu_representatives
from src.filtering import filter_by_min_genome_length
from src.taxonomy import extract_ictv_and_host_class, get_min_genome_lengths

# Load and process
df = load_metadata(config['data']['input_metadata'])
df = extract_ictv_and_host_class(df)

# Filter by class-specific minimum lengths
min_lengths = get_min_genome_lengths()
df_filtered = filter_by_min_genome_length(df, min_lengths)

# Analyze
print(f"Filtered to {df_filtered['votu'].nunique():,} vOTUs")
print("\nTop 10 viral classes:")
print(df_filtered['class'].value_counts().head(10))
```

## Output Files

### Tables
- `final_filtered_linear_complete.tsv` - Master dataset with all filters applied
- `linear_complete_high_quality.tsv` - High-quality vOTUs only
- `linear_complete_medium_quality.tsv` - Medium-quality vOTUs
- `summary_statistics.tsv` - Key metrics

### Figures
- `Fig2_overview.svg` - Dataset size comparison and distributions
- `Fig4_5_taxonomy.svg` - Taxonomy-based breakdown
- `Fig_SX_length_boxplots.svg` - Genome length by class

## Performance

Pipeline performance on full metaVR dataset (~12.7M records):
- Data loading: ~5 seconds
- Filtering & taxonomy: ~30 seconds
- Terminal analysis: ~2-3 minutes
- Visualization: ~1 minute

Total runtime: ~5-10 minutes (single core, 8GB RAM)

## Methods

### Filtering Strategy

1. **Topology segregation**: Identify exclusively linear vOTUs
2. **Nested region removal**: Exclude fragmented sequences in longer vOTUs
3. **Length filtering**: Keep UViGs within ±25bp of representative
4. **Group size**: Require minimum 3 UViGs per vOTU for statistical power
5. **Taxonomy filtering**: Remove sequences below class-specific length minimums

### Terminal Sequence Analysis

Terminal sequences (start/end 25bp) are compared within each vOTU to identify orthologs:
- Representative provides reference sequence
- Non-representative members compared for mismatches
- Allow up to 3 mismatches per terminus
- Require members from ≥3 unique samples

## Troubleshooting

**Issue: FileNotFoundError for data files**
- Verify paths in `config/config.yaml` are correct
- Use absolute paths, not relative paths

**Issue: Memory error with large datasets**
- Filter data before loading (use metadata pre-filtering)
- Process in chunks instead of all at once

**Issue: Conda/Mamba installation conflicts**
```bash
# Use conda-lock for reproducible environment
conda-lock create -f environment.yml
conda-lock install
```

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make changes with tests
4. Submit pull request

## Citation

If you use this pipeline, please cite:

```bibtex
@software{viral_complete_genomes_2024,
  title={Viral Linear Complete Genomes Analysis Pipeline},
  author={Your Name},
  url={https://github.com/yourusername/viral-complete-genomes},
  year={2024}
}
```

## License

This project is licensed under the MIT License - see LICENSE file for details.

## References

- [Viral Metadata Resource (metaVR)](https://zenodo.org/...)
- [ICTV Master Species List](https://ictv.global/)
- [Viral Sequence Classification](https://genomad.readthedocs.io/)

## Contact

Questions or issues? Open an issue on [GitHub](https://github.com/yourusername/viral-complete-genomes/issues)

---

**Last updated**: May 2024
**Version**: 1.0.0
