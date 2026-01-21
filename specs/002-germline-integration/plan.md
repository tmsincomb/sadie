# Implementation Plan: Germline Database Integration

**Branch**: `002-germline-integration` | **Date**: 2026-01-21 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/002-germline-integration/spec.md`

## Summary

Connect SADIE's new germline database module to existing `sadie.airr.Airr` and `sadie.renumbering.Renumbering` modules, enabling selection of germline providers (IMGT, OGRDB, VDJbase, custom) via a `germline_backend` parameter. The integration uses a feature flag (`SADIE_USE_GERMLINES_MODULE`) for backwards compatibility, with G3 API remaining the default while local germlines is opt-in. Test suite in `tests/unit/germlines/` mirrors critical path tests from airr and renumbering modules.

## Technical Context

**Language/Version**: Python 3.10+
**Primary Dependencies**: biopython (≥1.80), pyhmmer (^0.11.1), pydantic (≥2.0.0), pandas (≥1.5), PyYAML, click (≥8.0)
**Storage**: Local filesystem (FASTA files, BLAST databases, HMM models, YAML configs)
**Testing**: pytest with coverage
**Target Platform**: Linux/macOS (cross-platform Python)
**Project Type**: Single Python package with CLI
**Performance Goals**: Equivalent to G3 backend (NFR-001)
**Constraints**: No silent fallback to G3 on failure (NFR-002), offline-capable after initial setup
**Scale/Scope**: Supports human, mouse, rat, rabbit, macaque, dog species; V/D/J/C segments; H/K/L chains

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

The project constitution is a template without project-specific principles. Based on the codebase patterns observed:

| Principle | Status | Notes |
|-----------|--------|-------|
| Backwards Compatibility | ✅ PASS | Feature flag preserves G3 as default; germlines is opt-in |
| Test Coverage | ✅ PASS | Mirrored test suite in tests/unit/germlines/ |
| Local-First Operation | ✅ PASS | Germlines module eliminates G3 API dependency |
| Clear Error Messages | ✅ PASS | Fail with setup instructions, no silent fallback |

## Project Structure

### Documentation (this feature)

```text
specs/002-germline-integration/
├── spec.md              # Feature specification (complete)
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
├── contracts/           # Phase 1 output (Python API contracts)
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code (existing repository structure)

```text
src/sadie/
├── germlines/                    # NEW: Local germline database module
│   ├── __init__.py              # Public API (get_germline_genes, get_manager, etc.)
│   ├── models.py                # GermlineGene, ProviderMetadata (pydantic)
│   ├── manager.py               # GermlineManager (multi-source lookup)
│   ├── pipeline.py              # GermlinePipeline (3-stage processing)
│   ├── renumbering_integration.py  # LocalHMMBuilder
│   ├── g3_adapter.py            # Format conversion to G3 response format
│   ├── providers/               # Data source handlers
│   │   ├── base.py, imgt.py, ogrdb.py, vdjbase.py, custom.py
│   ├── builders/                # Data processing
│   │   ├── gapper.py, blast.py, aux.py, hmm.py
│   ├── utils/
│   │   └── feature_flags.py     # SADIE_USE_GERMLINES_MODULE flag
│   ├── sources/                 # Raw germline data (FASTA)
│   ├── normalized/              # Processed sequences
│   ├── igblast/                 # BLAST-ready databases
│   └── hmms/                    # Cached HMM models
│
├── airr/                         # EXISTING: AIRR annotation
│   ├── igblast/
│   │   └── germline.py          # ✅ INTEGRATED: GermlineData paths updated
│   └── ...
│
├── renumbering/                  # EXISTING: Antibody numbering
│   ├── aligners/
│   │   └── hmmer.py             # ✅ INTEGRATED: LocalHMMBuilder support
│   └── ...
│
└── reference/                    # EXISTING: Reference system
    └── reference.py             # ✅ INTEGRATED: GermlineManager + G3Adapter

tests/
├── unit/
│   ├── airr/                    # Existing AIRR tests (G3-based)
│   ├── renumbering/             # Existing renumbering tests (G3-based)
│   └── germlines/               # NEW: Mirrored tests using germlines backend
│       ├── test_airr_integration.py       # ✅ CREATED
│       ├── test_renumbering_integration.py # ✅ CREATED
│       └── test_reference_integration.py   # ✅ CREATED
└── integration/
    └── ...
```

**Structure Decision**: Existing single-package structure maintained. New `germlines/` module added under `src/sadie/`. Test suite mirrored in `tests/unit/germlines/`.

## Complexity Tracking

No constitution violations requiring justification. The integration follows established patterns:
- Feature flag for gradual migration (standard practice)
- Adapter pattern for format compatibility (G3Adapter)
- Provider pattern for multiple data sources (already implemented in germlines module)

## Integration Status (Current)

Based on codebase exploration, the integration is **67% complete** (34/51 tasks):

| Component | File | Status |
|-----------|------|--------|
| IgBLAST paths | `airr/igblast/germline.py` | ✅ Complete |
| Renumbering HMM | `renumbering/aligners/hmmer.py` | ✅ Complete |
| Reference system | `reference/reference.py` | ✅ Complete |
| Feature flag | `germlines/utils/feature_flags.py` | ✅ Complete |
| Test suite | `tests/unit/germlines/` | ✅ 82% Complete |
| BLAST databases | `germlines/igblast/` | ⏳ Blocked (needs BLAST+) |
| Documentation | INTEGRATION_GUIDE.md | ✅ Complete |

## Remaining Work

1. **BLAST database building** - Requires BLAST+ tools installation
2. **End-to-end validation** - Run full test suite with germlines backend
3. **Offline testing** - Verify no network calls after initial setup
4. **Result comparison** - Ensure germlines IMGT matches G3 IMGT output
