# Research: Germlines Module Completion

**Date**: 2026-01-14 | **Spec**: [spec.md](./spec.md) | **Plan**: [plan.md](./plan.md)

## Overview

This document captures research findings for the Phase 0 unknowns identified in the implementation plan. Each research question is addressed with findings, decisions, and implementation implications.

---

## 1. VDJbase Data Format & Access

**Status**: RESOLVED

### Question
What is VDJbase's current FASTA format and download mechanism?

### Findings

#### API Structure

**Base URL**: `https://vdjbase.org/admin/api/`

| Endpoint | Purpose | Example |
|----------|---------|---------|
| `genomic/species` | List supported species | Returns `["Rhesus Macaque", "Human"]` |
| `genomic/data_sets/{species}` | List datasets per species | IGH, IGK, IGL, IGHC for Human |
| `repseq/samples/{species}/{chain}` | Repertoire sample metadata | Paginated, includes genotype file paths |
| `repseq/sequences/{species}/{chain}` | Germline sequences | Returns allele sequences with metadata |

#### Supported Species & Datasets

**Human**:
- IGH (created 2025-07-09)
- IGK (created 2025-07-09)
- IGL (created 2025-07-08)
- IGHC - constant region (created 2025-08-06)

**Rhesus Macaque**:
- IGH (created 2025-09-28)
- IGK (created 2025-09-28)
- IGL (created 2025-09-28)

#### Sequence Data Format

The `/repseq/sequences/{species}/{chain}` endpoint returns paginated JSON:

```json
{
  "samples": [
    {
      "name": "IGHV2-70*15_g303a",
      "seq": "caggtcaccttgagggagtctggtcct...",
      "seq_len": "322",
      "type": "IGHV",
      "family": "IGHV2",
      "gene_name": "IGHV2-70",
      "species": "Human",
      "appears": 66,
      "novel": true,
      "is_single_allele": true,
      "low_confidence": true,
      "dataset": "IGH"
    }
  ]
}
```

**Key fields**:
- `name`: Allele identifier (mutation-based naming for novel alleles, e.g., `IGHV2-70*15_g303a`)
- `seq`: Ungapped nucleotide sequence
- `type`: Segment type (IGHV, IGHD, IGHJ)
- `novel`: Boolean - novel allele not in IMGT/OGRDB
- `low_confidence`: Quality indicator
- `appears`: Observation count across samples

#### OGRDB Report Files (Gapped Sequences)

Per-sample genotype files contain IMGT-gapped sequences:

**URL Pattern**:
```
https://vdjbase.org/static/study_data/VDJbase/samples/{species}/{chain}/{study}/{sample}/{sample}_ogrdb_report.csv
```

**CSV Format**:
```csv
"sequence_id","nt_sequence","nt_sequence_gapped"
"IGHV1-2*02","GCTTCTGGATACACCTTC...","NNNNNNN...GCTTCTGGATACACCTTC........."
```

#### Data Volume (Human IGH)

| Metric | Value |
|--------|-------|
| Total IGHV sequences | ~725 |
| Novel alleles | ~67% |
| Pagination | 250 items/page recommended |

#### API Quirks

1. **Page size affects results**: Smaller page sizes (10-100) recommended for consistent behavior
2. **Sorting varies**: Different page sizes may return different ordering
3. **No bulk FASTA download**: Must paginate through API

### Decision

**Approach**: Hybrid download with local caching

1. **Download Phase**: Paginate `/repseq/sequences` API, cache to local FASTA
2. **Gapped Templates**: Use OGRDB report CSVs or derive from IMGT gapped references
3. **Storage**: `src/sadie/germlines/sources/vdjbase/{species}/{chain}.fasta`
4. **Metadata**: Store `novel`, `low_confidence`, `appears` in FASTA headers

### Implementation Impact

```python
class VDJbaseProvider(GermlineProvider):
    """VDJbase integration using REST API with local caching."""
    
    BASE_URL = "https://vdjbase.org/admin/api"
    
    def fetch_genes(self, species: str, segment: str, chain: str) -> List[GermlineGene]:
        # 1. Check local cache first (offline-first)
        # 2. If not cached, paginate API and cache results
        # 3. Filter by segment type
        # 4. Note: sequences are UNGAPPED - require gapper module
        pass
```

---

## 2. Biopython Alignment Strategy for Gapping

**Status**: PARTIALLY RESOLVED

### Question
How to use Biopython alignment for auto-gapping sequences to IMGT format?

### Current Decision

Use per-gene IMGT-gapped templates when available; fallback to per-segment consensus template.

**Template source**: IMGT provider gapped FASTA in `src/sadie/germlines/sources/imgt/`

### Remaining Research Needed

1. **Identify available gapped references** in IMGT sources
2. **Evaluate Bio.Align/PairwiseAligner settings** for AA alignment:
   - Appropriate substitution matrix (BLOSUM62?)
   - Gap open/extend penalties
   - End gap handling
3. **Consensus template derivation** from gapped sequences

### Implementation Approach

```python
from Bio.Align import PairwiseAligner

class GapperService:
    def gap_sequence(self, sequence: str, segment: str, scheme: str = "imgt") -> str:
        # 1. Translate to amino acids
        # 2. Find best matching IMGT-gapped template
        # 3. Align AA sequences using PairwiseAligner
        # 4. Map gaps back to nucleotide positions
        # 5. Return gapped nucleotide sequence
        pass
```

---

## 3. G3 API Response Format

**Status**: NEEDS RESEARCH

### Question
What exact JSON format does G3 return for gene queries?

### Research Approach

1. Inspect existing G3 client code in `src/sadie/renumbering/clients/`
2. Capture sample API responses
3. Document field mappings for adapter pattern

### Expected Fields (from plan.md)

```json
{
  "gene": "string",
  "sequence": "string",
  "sequence_gapped": "string",
  "species": "string",
  "segment": "string",
  "chain": "string",
  "source": "string",
  "functional": "boolean",
  "regions": {
    "fwr1": "string",
    "cdr1": "string",
    "fwr2": "string",
    "cdr2": "string",
    "fwr3": "string"
  }
}
```

### Action Items

- [ ] Read `src/sadie/renumbering/clients/g3.py` to understand current implementation
- [ ] Document actual API response structure
- [ ] Identify all fields needed by Reference system

---

## 4. IgBLAST Auxiliary File Format

**Status**: NEEDS RESEARCH

### Question
What CDR/FWR annotations are required in auxiliary files?

### Research Approach

1. Review IgBLAST documentation
2. Inspect existing Sadie auxiliary files in `src/sadie/germlines/`
3. Document required columns and format

### Expected Format (preliminary)

TSV file with columns:
- Gene name
- FWR1 start/end
- CDR1 start/end
- FWR2 start/end
- CDR2 start/end
- FWR3 start/end

### Action Items

- [ ] Read existing `src/sadie/germlines/builders/aux.py`
- [ ] Review IgBLAST auxiliary file specifications
- [ ] Document format requirements

---

## 5. IMGT/OGRDB Download URLs

**Status**: NEEDS RESEARCH

### Question
Current stable URLs for bulk FASTA downloads from IMGT and OGRDB?

### IMGT (Preliminary)

- IMGT/GENE-DB: http://www.imgt.org/genedb/
- Bulk download may require programmatic parsing

### OGRDB

- API: https://ogrdb.airr-community.org/api/
- FASTA downloads available per species/locus

### Action Items

- [ ] Document stable IMGT download URLs
- [ ] Document OGRDB API endpoints for FASTA download
- [ ] Test download scripts for reliability
- [ ] Add retry/resume logic for large downloads

---

## Summary

| Research Question | Status | Priority |
|-------------------|--------|----------|
| VDJbase Data Format & Access | RESOLVED | - |
| Biopython Alignment Strategy | PARTIAL | High |
| G3 API Response Format | NEEDS RESEARCH | Medium |
| IgBLAST Auxiliary File Format | NEEDS RESEARCH | Medium |
| IMGT/OGRDB Download URLs | NEEDS RESEARCH | High |

## Next Steps

1. Complete remaining research items before implementation
2. Update plan.md with any changes to approach
3. Proceed to Phase 1 design artifacts once research complete
