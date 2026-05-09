# Getting Started Guide

## Quick Setup (5 minutes)

### 1. Clone and Install

```bash
# Clone the repository
git clone https://github.com/yourusername/viral-complete-genomes.git
cd viral-complete-genomes

# Create virtual environment
python -m venv venv
source venv/bin/activate          # Linux/Mac
# OR
venv\Scripts\activate             # Windows

# Install dependencies
pip install -r requirements.txt
```

### 2. Configure Data Paths

Edit `config/config.yaml`:

```yaml
data:
  input_metadata: "/path/to/your/UVIG_metadata.tsv"
  header_repr: "/path/to/representatives.txt"
  exclude_nested: "/path/to/excluded_sequences.txt"
  # Add other paths...
```

Replace `/path/to/` with actual paths to your data files.

### 3. Run Analysis

```bash
# Option A: Run complete pipeline
python scripts/run_pipeline.py --config config/config.yaml

# Option B: Run step-by-step examples
python examples/01_basic_filtering.py
python examples/02_terminal_sequence_analysis.py

# Option C: Use in your own Python code
python
>>> from src.filtering import filter_by_min_genome_length
>>> # ... your custom analysis
```

## Understanding the Workflow

```
Your Data (TSV)
      ↓
[1] Load metadata
      ↓
[2] Filter by topology (Linear only)
      ↓
[3] Remove nested sequences
      ↓
[4] Split by quality tier
      ↓
[5] Filter by length consistency (±25bp)
      ↓
[6] Require minimum group size (3+ UViGs)
      ↓
[7] Extract taxonomy (ICTV class)
      ↓
[8] Filter by class-specific minimum length
      ↓
[9] Analyze terminal sequences (optional)
      ↓
Final Dataset (TSV) + Summary Stats
```

## Common Tasks

### Task 1: Basic Filtering Only

```python
from src.data_loading import load_metadata, extract_topology_quality_subsets
from src.filtering import filter_exclusive_linear_votus

# Load
df = load_metadata('data.tsv')
subsets = extract_topology_quality_subsets(df)

# Get linear only
linear_df = subsets['all_linear']

# Filter
cleaned, n = filter_exclusive_linear_votus(linear_df)
print(f"{n} vOTUs have exclusively linear topology")
```

### Task 2: Get Specific Viral Classes

```python
from src.taxonomy import extract_ictv_and_host_class

df = extract_ictv_and_host_class(df)

# Filter to DNA viruses only
dna_viruses = df[df['class'].isin([
    'Caudoviricetes', 'Faserviricetes', 'Laserviricetes'
])]

print(f"Found {len(dna_viruses)} DNA virus sequences")
```

### Task 3: Analyze Host Associations

```python
from src.taxonomy import extract_taxonomy_level

# Extract host genus
df['host_genus'] = df['host_taxonomy'].apply(
    lambda x: extract_taxonomy_level(x, 'g__')
)

# Count by host
print(df['host_genus'].value_counts().head(10))
```

### Task 4: Find Genome Length Outliers

```python
from src.filtering import identify_length_outliers

outliers = identify_length_outliers(df, max_delta_bp=1000)
print(f"Found {len(outliers)} vOTUs with >1000bp length variation")
print(outliers[['votu', 'min_length', 'max_length', 'length_delta']].head())
```

## Data Requirements Checklist

✅ Metadata TSV with required columns:
- [ ] `uvig` - Unique genome identifier
- [ ] `votu` - Viral OTU cluster
- [ ] `length` - Genome length in bp
- [ ] `uvig_topology` - Genome topology
- [ ] `quality` - Quality assessment
- [ ] `ictv_taxonomy` - Viral taxonomy
- [ ] `host_taxonomy` - Host taxonomy (can be NaN)
- [ ] `source` - Data source

✅ Representative UViGs file:
- [ ] One UViG ID per line
- [ ] File format: text file

## Output Files Explained

After running the pipeline, check:

```
results/
├── tables/
│   ├── final_filtered_linear_complete.tsv     ← Main result
│   ├── linear_complete_high_quality.tsv       ← Quality-specific
│   ├── linear_complete_medium_quality.tsv
│   └── summary_statistics.tsv                 ← Key metrics
│
└── figures/                                    ← Visualizations
    ├── Fig2_overview.svg
    └── Fig4_5_taxonomy.svg
```

**Main result file** (`final_filtered_linear_complete.tsv`) contains:
- All original columns plus:
- `class` - Viral ICTV class
- `host_class` - Host taxonomy class
- `Genome_type` - DNA/RNA classification

## Customization Examples

### Change Filtering Strictness

In `config/config.yaml`:

```yaml
filtering:
  # Strict filtering (original)
  max_bp_difference: 25
  min_uvigs_per_votu: 3

  # OR Lenient filtering
  max_bp_difference: 100
  min_uvigs_per_votu: 2
```

Then re-run pipeline:
```bash
python scripts/run_pipeline.py --config config/config.yaml
```

### Focus on Specific Viral Classes

```python
from src.data_loading import load_metadata
from src.taxonomy import extract_ictv_and_host_class

df = load_metadata('data.tsv')
df = extract_ictv_and_host_class(df)

# Get only bacteriophages (DNA viruses)
phages = df[df['Genome_type'] == 'DNA']

# Get only Caudoviricetes (largest order of phages)
caudovirales = phages[phages['class'] == 'Caudoviricetes']

print(f"Found {len(caudovirales)} Caudoviricete sequences")
```

### Export Specific Subset

```python
# Filter to human-associated viruses
human_viruses = df[
    df['host_taxonomy'].str.contains('Homo sapiens', na=False)
]

# Save
human_viruses.to_csv('human_associated_viruses.tsv', sep='\t', index=False)
print(f"Exported {len(human_viruses)} human-associated sequences")
```

## Troubleshooting

### Problem: "FileNotFoundError: input_metadata"
**Solution**: Update paths in `config/config.yaml` to use absolute paths

```yaml
# ❌ Don't use relative paths
data:
  input_metadata: "data/metadata.tsv"

# ✅ Use absolute paths
data:
  input_metadata: "/Users/name/Documents/data/metadata.tsv"
```

### Problem: "ModuleNotFoundError: src"
**Solution**: Run from repository root directory

```bash
# ❌ Wrong
cd viral-complete-genomes/scripts
python run_pipeline.py

# ✅ Correct
cd viral-complete-genomes
python scripts/run_pipeline.py --config config/config.yaml
```

### Problem: Memory error with large dataset
**Solution**: Process in chunks or filter before loading

```python
# Method 1: Chunked processing
chunksize = 100000
for chunk in pd.read_csv('data.tsv', sep='\t', chunksize=chunksize):
    # Process chunk
    filtered = process(chunk)
    results.append(filtered)

# Method 2: Pre-filter with shell tools
# grep "Linear" data.tsv > linear_only.tsv
# Then load linear_only.tsv
```

## Next Steps

1. **Read the STRUCTURE.md** for detailed module documentation
2. **Run examples/** for hands-on tutorials
3. **Adapt config/config.yaml** for your data
4. **Execute scripts/run_pipeline.py** for full analysis
5. **Check results/** for outputs

## Getting Help

- Check **README.md** for overview
- Check **STRUCTURE.md** for code documentation
- Review **examples/** for working code samples
- See docstrings: `help(function_name)` in Python
- Open GitHub issue with error details

## Next: Adapt for Your Data

### Step 1: Prepare your metadata
Ensure your TSV has required columns listed in "Data Requirements Checklist"

### Step 2: Update config
Edit `config/config.yaml` with your file paths

### Step 3: Run example
```bash
python examples/01_basic_filtering.py
```

### Step 4: Verify output
Check `results/example_filtered_high_quality.tsv`

### Step 5: Run full pipeline
```bash
python scripts/run_pipeline.py --config config/config.yaml
```

---

**You're ready to begin!** 🚀

Questions? Open an issue on GitHub or check README.md
