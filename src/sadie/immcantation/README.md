# SADIE Immcantation Pipeline

## Overview

The SADIE Immcantation Pipeline provides a comprehensive framework for analyzing VDJ antibody repertoires from Illumina sequencing data. This module integrates with the Immcantation ecosystem (pRESTO, Change-O, Dowser) while leveraging SADIE's native capabilities for robust, end-to-end analysis.

## Features

### Complete Workflow Support

1. **pRESTO Integration** - Quality control and preprocessing
   - Quality and length filtering
   - Primer identification and masking
   - UMI consensus building
   - Paired-end assembly
   - Deduplication

2. **Change-O Integration** - VDJ assignment and clonal clustering
   - IgBLAST-based V(D)J gene assignment
   - AIRR-compliant output
   - Distance-based clonal clustering
   - Germline sequence reconstruction

3. **Lineage Analysis** - Phylogenetic tree construction
   - Clone size filtering
   - Neighbor-joining tree construction
   - Tree statistics and metrics
   - Support for multiple tree building methods

4. **SADIE Native Analysis** - Enhanced capabilities
   - Antibody numbering schemes (IMGT, Kabat, Chothia, Martin, Aho)
   - Mutational analysis
   - Inferred germline assignment
   - CDR-based clustering

## Installation

### Prerequisites

The Immcantation pipeline requires:
- Python 3.9+
- SADIE package
- Immcantation tools (pRESTO, Change-O) - optional
- Biopython for tree building

### Quick Install

```bash
# Install SADIE
pip install -e /path/to/sadie

# Optional: Install Immcantation tools
pip install presto changeo
```

## Usage

### Example 1: Downstream Analysis (Pre-annotated Data)

For data that's already been annotated with VDJ genes:

```python
from pathlib import Path
import pandas as pd
from sadie.airr.airrtable import AirrTable
from sadie.cluster import Cluster
from sadie.immcantation.lineage import LineageAnalyzer

# Load pre-annotated data
airr_table = pd.read_feather("annotated_data.feather")
airr_table = AirrTable(airr_table)

# Filter productive sequences
productive = airr_table[airr_table["productive"] == True]

# Clonal clustering
cluster_api = Cluster(
    productive,
    linkage="single",
    lookup=["cdr1_aa", "cdr2_aa", "cdr3_aa"],
)
clustered = cluster_api.cluster(0.15)  # 15% distance threshold

# Build lineage trees
lineage = LineageAnalyzer()
filtered_clones = lineage.filter_clones_by_size(
    clustered, min_size=5, max_clone_size=100
)
clone_files = lineage.split_clones(filtered_clones, Path("./clones"))
trees = lineage.build_trees_simple(clone_files, Path("./trees"))
```

### Example 2: Complete Pipeline

For raw FASTA data:

```python
from sadie.immcantation.pipeline import ImmcantationPipeline
from sadie.immcantation.config import PipelineConfig

# Configure pipeline
config = PipelineConfig(
    output_dir=Path("./results"),
    run_presto=False,  # Skip if data is pre-processed
    run_changeo=True,
    run_lineage=True,
)

# Run pipeline
pipeline = ImmcantationPipeline(config)
outputs = pipeline.run(
    input_file=Path("sequences.fasta"),
    sample_name="my_sample",
)
```

### Example 3: Using SADIE Native Tools

```python
from sadie.airr import Airr
from sadie.airr.methods import run_mutational_analysis
from sadie.cluster import Cluster

# VDJ annotation
airr_api = Airr("human")
airr_table = airr_api.run_fasta("sequences.fasta")

# Mutational analysis
airr_table = run_mutational_analysis(airr_table, scheme="imgt")

# Clustering
cluster_api = Cluster(airr_table, lookup=["cdr3_aa"])
clustered = cluster_api.cluster(0.10)
```

## Command-Line Usage

### Downstream Analysis Demo

```bash
python examples/downstream_analysis_demo.py
```

This demonstrates:
- Loading pre-annotated AIRR data
- Clonal clustering
- Clone statistics
- Lineage tree construction
- Output file generation

### Full Pipeline

```bash
python examples/run_immcantation_pipeline.py \
    --input sequences.fasta \
    --output results/ \
    --organism human \
    --distance-threshold 0.15 \
    --min-clone-size 5
```

## Configuration

### Pipeline Configuration

```python
from sadie.immcantation.config import (
    PipelineConfig,
    PrestoConfig,
    ChangeoConfig,
    LineageConfig,
)

# Detailed configuration
config = PipelineConfig(
    presto=PrestoConfig(
        min_quality=20,
        min_length=200,
        nproc=4,
    ),
    changeo=ChangeoConfig(
        organism="human",
        distance_threshold=0.15,
        linkage_method="single",
        nproc=4,
    ),
    lineage=LineageConfig(
        min_clone_size=5,
        max_clone_size=1000,
        method="igphyml",
    ),
    output_dir=Path("./results"),
    verbose=True,
)
```

### Primer Sets

Pre-configured primer sets for common protocols:

```python
from sadie.immcantation.config import PRIMER_SETS

# Human heavy chain 5'RACE primers
heavy_primers = PRIMER_SETS["human_heavy_5race"]

# Kappa light chain primers
kappa_primers = PRIMER_SETS["human_kappa_5race"]

# Lambda light chain primers
lambda_primers = PRIMER_SETS["human_lambda_5race"]
```

## Output Files

### Standard Outputs

1. **Annotation** (`*_db-pass.tsv`)
   - AIRR-compliant V(D)J annotations
   - Productive sequence flags
   - CDR/FWR sequences

2. **Clonal Clustering** (`*_clone-pass.tsv`)
   - Clone ID assignments
   - Distance-based clustering results

3. **Clone Statistics** (`*_clone_stats.tsv`)
   - Clone sizes
   - V/J gene usage per clone
   - Junction length statistics

4. **Lineage Trees** (`trees/*.nwk`)
   - Newick format phylogenetic trees
   - One tree per clonal family

5. **Reports** (`*_report.txt`)
   - Pipeline summary
   - Statistics
   - Output file listing

## Architecture

### Module Structure

```
sadie/immcantation/
├── __init__.py           # Package initialization
├── config.py             # Configuration classes
├── pipeline.py           # Main pipeline orchestrator
├── presto_wrapper.py     # pRESTO command wrappers
├── changeo_wrapper.py    # Change-O command wrappers
├── lineage.py            # Lineage analysis tools
└── README.md             # This file
```

### Pipeline Flow

```
Input FASTA/FASTQ
       ↓
   [pRESTO]              (Optional)
   - Quality Control
   - UMI Processing
   - Consensus Building
       ↓
   [Change-O / SADIE]
   - V(D)J Assignment
   - Productive Filtering
       ↓
   [Clustering]
   - CDR-based Distance
   - Hierarchical Clustering
       ↓
   [Lineage Analysis]
   - Clone Filtering
   - Tree Construction
   - Statistics
       ↓
   Output Files
```

## Performance

### Optimization Tips

1. **Use Multi-Processing**
   ```python
   config.changeo.nproc = 8
   config.presto.nproc = 8
   ```

2. **Filter Early**
   ```python
   # Filter by quality/length before annotation
   productive_only = airr_table[airr_table["productive"] == True]
   ```

3. **Adjust Clone Size Limits**
   ```python
   # Skip very large clones for tree building
   config.lineage.max_clone_size = 100
   ```

### Benchmarks

On a modern 8-core CPU:
- **Annotation**: ~1000 sequences/minute (SADIE IgBLAST)
- **Clustering**: ~10000 sequences/minute
- **Tree Building**: ~50 trees/minute (NJ, clone size < 50)

## Troubleshooting

### Common Issues

**IgBLAST not found**
```python
# Specify path explicitly
airr_api = Airr("human", igblast_exe="/path/to/igblastn")
```

**Change-O command not found**
```bash
# Activate virtual environment first
source /path/to/venv/bin/activate
python examples/run_pipeline.py ...
```

**Empty clone trees**
```python
# Lower minimum clone size
config.lineage.min_clone_size = 3
```

**Memory issues with large datasets**
```python
# Process in batches
for batch in batch_sequences(large_file, batch_size=1000):
    process_batch(batch)
```

## Citation

If you use this pipeline in your research, please cite:

- **SADIE**: (SADIE publication)
- **pRESTO**: Vander Heiden JA, et al. (2014) Bioinformatics 30(13):1930-1932
- **Change-O**: Gupta NT, et al. (2015) Bioinformatics 31(20):3356-3358
- **IgBLAST**: Ye J, et al. (2013) Nucleic Acids Res 41:W34-W40

## Support

- GitHub Issues: https://github.com/jwillis0720/sadie/issues
- Documentation: https://sadie.readthedocs.io
- Immcantation Docs: https://immcantation.readthedocs.io

## License

MIT License - see LICENSE file for details.
