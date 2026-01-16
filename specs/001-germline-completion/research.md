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

**Status**: RESOLVED

### Question
How to use Biopython alignment for auto-gapping sequences to IMGT format?

### Decision

Use per-gene IMGT-gapped templates when available; fallback to per-segment consensus template.

**Template source**: IMGT provider gapped FASTA in `src/sadie/germlines/sources/imgt/`

### Implementation (Completed)

The `GapperService` has been implemented in `src/sadie/germlines/builders/gapper.py` with:

1. **BioPython PairwiseAligner** configured for antibody sequences:
   - Global alignment mode
   - Match score: 2.0, mismatch: -1.0
   - Gap open: -10.0, extend: -0.5
   - End gaps allowed without penalty

2. **Template Strategy**:
   - Per-gene template lookup (exact match by gene name)
   - Per-segment consensus fallback (built from available gapped sequences)
   - Templates cached for performance

3. **Workflow**:
   - Translate nucleotide to amino acid
   - Align AA sequences using PairwiseAligner
   - Extract gap positions from template
   - Apply gaps to nucleotide sequence at codon boundaries
   - Use "." (period) as gap character per IMGT convention

4. **Error Handling**:
   - D segments return None (no gapping required)
   - Failed gapping logs WARNING and returns None
   - Graceful degradation on missing templates

---

## 3. G3 API Response Format

**Status**: RESOLVED

### Question
What exact JSON format does G3 return for gene queries?

### Findings

From `src/sadie/renumbering/clients/g3.py`, the G3 client returns gene data with the following structure:

```json
{
  "gene": "IGHV1-69*01",
  "sequence": "CAGGTGCAGCTGGTGGAG...",
  "sequence_gapped": "CAGGTGCAGCTGGTGGAG......",
  "species": "human",
  "segment": "V",
  "chain": "H",
  "source": "imgt",
  "functional": true,
  "regions": {
    "fwr1": "CAGGTGCAGCTGGTGGAG",
    "cdr1": "GGTGGCAGCTTC",
    "fwr2": "TGGGTGCGCCAG",
    "cdr2": "ATAGACAGCAGTGGC",
    "fwr3": "CGCTCCGTGAAGGGCCGATTC"
  }
}
```

### Key Implementation Details

- G3 API endpoint: `https://g3.jordanrwillis.com/api/v1/genes`
- Uses LRU caching for performance
- Filters out certain species (pig, cow, cat, alpaca, rhesus, dog)
- Builds Stockholm alignments for HMM creation via `pyhmmer`

### Adapter Pattern

The Reference system adapter (`FR-012b`) should:
1. Call `GermlineManager.get_gene_by_name(name, species)`
2. Transform `GermlineGene` object to G3 response format
3. Return dict matching the structure above

---

## 4. IgBLAST Auxiliary File Format

**Status**: RESOLVED

### Question
What CDR/FWR annotations are required in auxiliary files?

### Findings

From `src/sadie/germlines/builders/aux.py`, the auxiliary file format is:

**Format**: Tab-separated values (TSV)

**Columns**:
1. `gene_name` - Gene identifier (e.g., "IGHV1-69*01")
2. `chain` - Chain type (H, K, L)
3. `segment` - Segment type (V, J)
4. `fwr1_start`, `fwr1_end` - Framework 1 boundaries
5. `cdr1_start`, `cdr1_end` - CDR1 boundaries
6. `fwr2_start`, `fwr2_end` - Framework 2 boundaries
7. `cdr2_start`, `cdr2_end` - CDR2 boundaries
8. `fwr3_start`, `fwr3_end` - Framework 3 boundaries

**IMGT Position Numbering** (from spec FR-037b):
- FWR1: positions 1-26
- CDR1: positions 27-38
- FWR2: positions 39-55
- CDR2: positions 56-65
- FWR3: positions 66-104

### Current Implementation Status

The `AuxFileBuilder` in `aux.py` is currently a stub. Key methods that need completion:
- `_parse_imgt_regions()` - Convert IMGT gaps to region boundaries
- `_create_aux_entry()` - Generate actual aux file entries

Output location: `igblast/aux_db/{scheme}/{species}_gl.aux`

---

## 5. IMGT/OGRDB Download URLs

**Status**: PARTIALLY RESOLVED

### Question
Current stable URLs for bulk FASTA downloads from IMGT and OGRDB?

### IMGT

- **IMGT/GENE-DB**: http://www.imgt.org/genedb/
- **Reference Directory**: http://www.imgt.org/download/GENE-DB/
- Bulk download requires programmatic parsing of HTML or FTP access
- Current `download_imgt.py` is a stub with manual instructions

### OGRDB

- **API Base**: https://ogrdb.airr-community.org/api/
- **API Documentation**: https://ogrdb.airr-community.org/api/docs
- FASTA downloads available per species/locus via API endpoints
- Current `OGRDBProvider` reads pre-downloaded files

### VDJbase (from Section 1)

- **API Base**: https://vdjbase.org/admin/api/
- **Sequences Endpoint**: `/repseq/sequences/{species}/{chain}`
- Returns paginated JSON with sequence data
- Manual FASTA download recommended (API pagination quirks)

### Remaining Action Items

- [ ] Implement automated IMGT download script
- [ ] Implement automated OGRDB download script
- [ ] Add retry/resume logic for large downloads
- [ ] Test download scripts for reliability

---

## Summary

| Research Question | Status | Priority |
|-------------------|--------|----------|
| VDJbase Data Format & Access | RESOLVED | - |
| Biopython Alignment Strategy | RESOLVED | - |
| G3 API Response Format | RESOLVED | - |
| IgBLAST Auxiliary File Format | RESOLVED | - |
| IMGT/OGRDB Download URLs | PARTIALLY RESOLVED | Medium |

## Next Steps

1. ~~Complete remaining research items before implementation~~ **DONE** (mostly)
2. Implement download scripts for IMGT and OGRDB
3. Complete `AuxFileBuilder` stub with region parsing
4. Proceed with implementation per tasks.md
