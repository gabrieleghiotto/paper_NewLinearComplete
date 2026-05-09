# Repository Structure

```
viral-complete-genomes/
│
├── README.md                          # Project overview and installation
├── STRUCTURE.md                       # This file - directory organization
├── requirements.txt                   # Python dependencies
├── LICENSE                            # MIT License
│
├── config/
│   └── config.yaml                    # Configuration file (EDIT THIS)
│                                      # - Data paths
│                                      # - Filter parameters
│                                      # - Taxonomy settings
│
├── src/                               # Main source code modules
│   ├── __init__.py
│   ├── data_loading.py               # Load and organize metadata
│   │   ├── load_metadata()
│   │   ├── load_representative_uvigs()
│   │   ├── extract_topology_quality_subsets()
│   │   ├── extract_quality_tiers()
│   │   ├── get_votu_representatives()
│   │   └── save_results()
│   │
│   ├── filtering.py                  # Filtering functions
│   │   ├── filter_exclusive_linear_votus()
│   │   ├── remove_nested_votus()
│   │   ├── filter_by_length_difference()
│   │   ├── filter_by_group_size()
│   │   ├── filter_linear_uvigs_comprehensive()
│   │   ├── filter_by_min_genome_length()
│   │   ├── identify_length_outliers()
│   │   └── print_filter_summary()
│   │
│   ├── taxonomy.py                   # Taxonomy classification
│   │   ├── extract_taxonomy_level()
│   │   ├── extract_ictv_and_host_class()
│   │   ├── extract_order_family()
│   │   ├── classify_genome_type()
│   │   ├── add_genome_type_column()
│   │   ├── get_min_genome_lengths()
│   │   └── get_class_summary()
│   │
│   └── terminal_analysis.py          # Terminal sequence analysis
│       ├── flag_representative_uvigs()
│       ├── compute_sequence_mismatches()
│       ├── compute_uvig_mismatches_per_member()
│       ├── filter_by_passing_members()
│       └── summarize_passing_uvigs()
│
├── scripts/
│   └── run_pipeline.py               # Complete analysis pipeline
│       Usage: python scripts/run_pipeline.py --config config/config.yaml
│
├── examples/
│   ├── __init__.py
│   ├── 01_basic_filtering.py        # Simple filtering workflow
│   │   - Load data
│   │   - Extract subsets
│   │   - Filter exclusively linear
│   │   - Apply quality filters
│   │   - Summary statistics
│   │
│   └── 02_terminal_sequence_analysis.py  # Terminal sequence analysis
│       - Load conserved sequences
│       - Compute mismatches
│       - Filter by conservation
│       - Generate statistics
│
├── results/                          # Output directory (created by scripts)
│   ├── tables/
│   │   ├── final_filtered_linear_complete.tsv
│   │   ├── linear_complete_high_quality.tsv
│   │   └── summary_statistics.tsv
│   ├── figures/
│   │   ├── Fig2_overview.svg
│   │   └── Fig4_5_taxonomy.svg
│   └── filtered_uvigs/
│       ├── uvig_linear_complete.txt
│       └── uvig_new_linear_complete.txt
│
└── notebooks/                        # Jupyter notebook examples (optional)
    ├── 01_data_loading.ipynb
    ├── 02_filtering.ipynb
    ├── 03_terminal_analysis.ipynb
    └── 04_visualization.ipynb
```

## Module Quick Reference

### data_loading.py
**Purpose**: Load and organize viral metadata

Key functions:
- `load_metadata()` - TSV → DataFrame
- `extract_topology_quality_subsets()` - Split by topology & quality
- `extract_quality_tiers()` - Extract High/Medium/ND tiers
- `save_results()` - DataFrame → TSV

### filtering.py
**Purpose**: Filter sequences by various criteria

Key functions:
- `filter_exclusive_linear_votus()` - Keep only linear vOTUs
- `filter_by_length_difference()` - ±bp from representative
- `filter_by_group_size()` - Min UViGs per vOTU
- `filter_linear_uvigs_comprehensive()` - Combined pipeline
- `filter_by_min_genome_length()` - Class-specific minimums

### taxonomy.py
**Purpose**: Extract and classify viral taxonomy

Key functions:
- `extract_ictv_and_host_class()` - Parse class from taxonomy
- `extract_order_family()` - Extract order/family ranks
- `classify_genome_type()` - DNA vs RNA classification
- `get_min_genome_lengths()` - Return class length dict

### terminal_analysis.py
**Purpose**: Analyze conserved terminal sequences

Key functions:
- `compute_uvig_mismatches_per_member()` - Mismatch counting
- `filter_by_passing_members()` - Filter by conservation
- `summarize_passing_uvigs()` - Generate statistics

## Typical Workflow

### Minimal Pipeline (5 minutes)
```
config/config.yaml (update paths)
         ↓
scripts/run_pipeline.py
         ↓
results/tables/final_filtered_linear_complete.tsv
```

### Step-by-Step Analysis
```
examples/01_basic_filtering.py
         ↓
examples/02_terminal_sequence_analysis.py
         ↓
Custom visualization notebooks
```

### Custom Analysis
```python
# Import specific modules as needed
from src.filtering import filter_by_min_genome_length
from src.taxonomy import extract_ictv_and_host_class

df = pd.read_csv('data.tsv', sep='\t')
df = extract_ictv_and_host_class(df)
df = filter_by_min_genome_length(df, min_len_dict)
```

## Configuration File Hierarchy

```
config/config.yaml
├── data: paths to input files
├── output: paths for results
├── filtering: parameter values
│   ├── max_bp_difference: 25
│   ├── min_uvigs_per_votu: 3
│   └── min_unique_samples: 3
├── terminal_analysis: mismatch thresholds
├── taxonomy: class-specific minimums
└── visualization: figure settings
```

## File Format Reference

### Input Metadata (TSV)
Required columns:
- `uvig`: Genome ID
- `votu`: Viral OTU ID
- `length`: Genome length (bp)
- `uvig_topology`: Linear/DTR/ITR
- `quality`: High/Medium/Not-determined
- `ictv_taxonomy`: Semicolon-delimited ICTV string
- `host_taxonomy`: Semicolon-delimited host string
- `source`: Metagenome/Metatranscriptome/Isolate

Optional columns:
- `genomad_score`
- `viral_confidence`
- `n_virus_hallmarks`

### Terminal Sequences (TSV)
Required columns:
- `full_name`: Full sequence identifier
- `start_25bp`: First 25bp
- `end_25bp`: Last 25bp

### Output Tables (TSV)
Standard columns from input plus:
- `class`: Extracted ICTV class
- `host_class`: Extracted host class
- `Genome_type`: DNA/RNA classification
- `order`: ICTV order (if extracted)
- `family`: ICTV family (if extracted)

## Performance Notes

- **Data loading**: ~5 sec (full metaVR: 12.7M records)
- **Filtering**: ~30 sec
- **Terminal analysis**: ~2-3 min
- **Figures**: ~1 min
- **Total**: ~5-10 min (single core, 8GB RAM)

## Extending the Pipeline

### Adding a New Filter
1. Create function in `src/filtering.py`:
```python
def my_new_filter(df, parameter):
    return df[df.column > parameter].copy()
```

2. Update `scripts/run_pipeline.py`:
```python
from src.filtering import my_new_filter
df = my_new_filter(df, config['filtering']['my_parameter'])
```

3. Test with examples:
```python
python examples/01_basic_filtering.py
```

### Adding Visualizations
1. Create function in `src/visualization.py` (new module):
```python
def plot_metric(df, output_path):
    fig, ax = plt.subplots()
    # ... plotting code ...
    fig.savefig(output_path)
```

2. Call from pipeline:
```python
from src.visualization import plot_metric
plot_metric(df, output_path)
```

## Troubleshooting

**File not found errors**
- Check paths in `config/config.yaml`
- Use absolute paths, not relative

**Memory errors**
- Process data in chunks
- Pre-filter data before loading

**Incorrect results**
- Verify input data format matches specification
- Check filter parameters in config

---

Last updated: May 2024
For questions, see README.md or open an issue on GitHub
