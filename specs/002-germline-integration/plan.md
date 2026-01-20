# Implementation Plan: Germline Database Integration

**Branch**: `002-germline-integration` | **Date**: 2026-01-20 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `/specs/002-germline-integration/spec.md`

**Note**: This template is filled in by the `/speckit.plan` command. See `.specify/templates/commands/plan.md` for the execution workflow.

## Summary

Connect SADIE's germline database module to existing AIRR annotation and renumbering systems. Enable users to select germline providers (IMGT, OGRDB, VDJbase, custom) instead of being limited to G3 API. Create mirrored test suite validating germlines backend produces equivalent results. Maintain backwards compatibility with G3 as default.

**Key Technical Approach**: Feature-flag controlled integration at three points (IgBLAST paths, HMM generation, Reference system) with adapter pattern for G3 API compatibility.

## Technical Context

**Language/Version**: Python 3.10+
**Primary Dependencies**: pyhmmer (HMM building), Biopython (sequence parsing), pydantic (models), existing sadie modules
**Storage**: File-based (FASTA, BLAST databases, HMM binaries, Stockholm alignments)
**Testing**: pytest with existing fixtures from airr/renumbering test suites
**Target Platform**: Linux/macOS (same as existing SADIE)
**Project Type**: Single (existing Python package integration)
**Performance Goals**: Equivalent to G3 API (<2s for gene lookup, <10s for HMM generation/caching)
**Constraints**: Zero breaking changes to existing APIs, offline-capable after initial setup, graceful fallback to G3
**Scale/Scope**: 3 integration points, ~500 lines new code, mirror ~20 critical path tests

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### ✅ Principle V: Integration Compatibility (PRIMARY GATE)

**Requirements**:
- Must maintain backward compatibility with existing Sadie components
- Migration path required for breaking changes

**Status**: ✅ PASS
- All integrations use feature flags (SADIE_USE_GERMLINES_MODULE env var)
- G3 API remains default; germlines is opt-in
- Graceful fallback on errors
- No breaking changes to public APIs

### ✅ Principle III: Local-First Operation

**Requirements**:
- All data stored locally after initial download
- No runtime dependency on external APIs

**Status**: ✅ PASS
- Germlines module already implements local-first (from 001-germline-completion)
- Integration leverages existing local data
- G3 fallback only for backwards compat, not required

### ✅ Principle I: Provider-Based Architecture

**Requirements**:
- Use provider pattern for data sources
- Providers self-contained and stateless

**Status**: ✅ PASS
- Integration uses existing GermlineManager with provider system
- No new providers needed (IMGT, OGRDB, VDJbase, custom already implemented)
- Stateless adapter pattern for G3 compatibility

## Project Structure

### Documentation (this feature)

```text
specs/002-germline-integration/
├── plan.md              # This file (/speckit.plan command output)
├── research.md          # Phase 0 output (/speckit.plan command)
├── data-model.md        # Phase 1 output (/speckit.plan command)
├── quickstart.md        # Phase 1 output (/speckit.plan command)
├── contracts/           # Phase 1 output (/speckit.plan command)
└── tasks.md             # Phase 2 output (/speckit.tasks command - NOT created by /speckit.plan)
```

### Source Code (repository root)

```text
src/sadie/
├── germlines/                      # Germlines module (from 001-germline-completion)
│   ├── __init__.py                # Public API
│   ├── models.py                  # GermlineGene, ProviderMetadata
│   ├── manager.py                 # GermlineManager
│   ├── pipeline.py                # Processing pipeline
│   ├── providers/                 # IMGT, OGRDB, VDJbase, custom
│   ├── builders/                  # BLAST, HMM, auxiliary file builders
│   ├── g3_adapter.py             # NEW: GermlineGene → G3 format adapter
│   └── renumbering_integration.py # NEW: LocalHMMBuilder for renumbering
│
├── airr/
│   └── igblast/
│       └── germline.py            # MODIFY: Add feature flag + path switching
│
├── reference/
│   └── reference.py               # MODIFY: Add germlines backend option
│
└── renumbering/
    └── aligners/
        └── hmmer.py               # MODIFY: Add LocalHMMBuilder integration

tests/
├── unit/
│   ├── airr/                      # Existing AIRR tests (unchanged)
│   ├── renumbering/               # Existing renumbering tests (unchanged)
│   └── germlines/                 # NEW: Mirrored integration tests
│       ├── test_airr_integration.py
│       └── test_renumbering_integration.py
```

**Structure Decision**: Single project integration (Option 1). Adding integration layer to existing SADIE package structure. Three integration points at IgBLAST paths, HMM generation, and Reference system. New germlines test directory mirrors critical paths from airr/renumbering tests.

## Complexity Tracking

No constitution violations - all gates pass with current design.

