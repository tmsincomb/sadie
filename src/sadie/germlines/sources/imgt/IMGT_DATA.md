# IMGT Germline Data

This directory should contain IMGT germline sequences downloaded from https://www.imgt.org/

## Directory Structure

```
imgt/
├── human/
│   ├── IGHV.fasta
│   ├── IGHD.fasta
│   ├── IGHJ.fasta
│   ├── IGKV.fasta
│   ├── IGKJ.fasta
│   ├── IGLV.fasta
│   └── IGLJ.fasta
├── mouse/
│   └── ... (same structure)
└── ... (other species)
```

## How to Obtain IMGT Data

### Method 1: Manual Download from IMGT

1. Visit https://www.imgt.org/vquest/refseqh.html
2. Select species (e.g., "Homo sapiens")
3. Download reference sequences for each segment:
   - Heavy chain V: IG Heavy V genes
   - Heavy chain D: IG Heavy D genes
   - Heavy chain J: IG Heavy J genes
   - Kappa chain V: IG Kappa V genes
   - Kappa chain J: IG Kappa J genes
   - Lambda chain V: IG Lambda V genes
   - Lambda chain J: IG Lambda J genes

4. Save each as FASTA format
5. Rename to match naming convention: `IG{chain}{segment}.fasta`
   - Example: Heavy V → `IGHV.fasta`
   - Example: Kappa J → `IGKJ.fasta`

6. Place in appropriate species directory

### Method 2: Automated Download (TODO)

```bash
# Future: automated download script
python scripts/download_imgt.py --species human mouse
```

## FASTA Format Expected

IMGT FASTA files have this format:

```
>IGHV1-2*01|Homo sapiens|F|...
cag.gtgcagctggtgcag...tctggggctgag...gtgaag...
```

**Key Features:**
- Headers contain gene name, species, functionality code
- Sequences are **IMGT-gapped** (dots indicate gaps)
- Sequences follow IMGT numbering scheme

## Important Notes

- **Gapped sequences**: IMGT provides sequences with gaps (dots)
- **Functionality**: "F" = functional, "ORF" = open reading frame, "P" = pseudogene
- **Multiple alleles**: Each gene may have multiple alleles (e.g., *01, *02)

## Validation

After adding files, validate with:

```python
from sadie.germlines import get_manager

manager = get_manager()
genes = manager.get_genes("human", "V", "H")
print(f"Loaded {len(genes)} IGHV genes from IMGT")
```

## References

- IMGT Homepage: https://www.imgt.org/
- IMGT/GENE-DB: https://www.imgt.org/genedb/
- Reference Sequences: https://www.imgt.org/vquest/refseqh.html
