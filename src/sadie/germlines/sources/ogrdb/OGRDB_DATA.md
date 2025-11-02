# OGRDB Germline Data

This directory should contain OGRDB germline sequences from https://ogrdb.airr-community.org/

## Directory Structure

```
ogrdb/
├── human/
│   ├── IGHV.fasta
│   ├── IGHD.fasta (if available)
│   ├── IGHJ.fasta
│   ├── IGKV.fasta
│   ├── IGKJ.fasta
│   ├── IGLV.fasta
│   └── IGLJ.fasta
└── ... (other species as available)
```

## How to Obtain OGRDB Data

### Method 1: Manual Download from OGRDB

1. Visit https://ogrdb.airr-community.org/
2. Browse germline sets
3. Select species and locus
4. Download sequences in FASTA format
5. Save with naming convention: `IG{chain}{segment}.fasta`

### Method 2: OGRDB API (TODO)

```python
# Future: automated download via API
import requests

base_url = "https://ogrdb.airr-community.org/api"
response = requests.get(f"{base_url}/germline/species/human/IGHV")
# Process and save...
```

### Method 3: Automated Script (TODO)

```bash
# Future: automated download script
python scripts/download_ogrdb.py --species human
```

## FASTA Format Expected

OGRDB FASTA format:

```
>IGHV1-69*01
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCT...
```

**Key Features:**
- Headers contain gene name
- Sequences may be **ungapped** or **gapped**
- Novel alleles discovered through repertoire sequencing

## Priority with IMGT

By default, OGRDB has lower priority than IMGT:
- If a gene exists in both → IMGT version is used
- Novel genes only in OGRDB → included
- Priority order: `custom > imgt > ogrdb`

To prioritize OGRDB over IMGT:

```python
from sadie.germlines import GermlineManager

manager = GermlineManager(providers=["custom", "ogrdb", "imgt"])
genes = manager.get_genes("human", "V", "H")
```

## Important Notes

- **Community-curated**: OGRDB sequences are community-submitted
- **Novel alleles**: May contain alleles not in IMGT
- **Evidence-based**: Each sequence has supporting evidence from repertoire data
- **Ungapped sequences**: OGRDB sequences may need gapping (automatic)

## Validation

After adding files:

```python
from sadie.germlines import get_manager

manager = get_manager()
genes = manager.get_genes("human", "V", "H")

# Check how many came from OGRDB
ogrdb_genes = [g for g in genes if g.source == "ogrdb"]
print(f"Loaded {len(ogrdb_genes)} genes from OGRDB")
```

## References

- OGRDB Homepage: https://ogrdb.airr-community.org/
- OGRDB API Docs: https://ogrdb.airr-community.org/api/docs
- AIRR Community: https://www.antibodysociety.org/the-airr-community/
- Publication: https://academic.oup.com/nar/article/48/D1/D964/5576123
