# Implementation Tasks: Germline Database Integration

**Feature**: 002-germline-integration
**Branch**: `002-germline-integration`
**Spec**: [spec.md](./spec.md) | **Plan**: [plan.md](./plan.md)

## Summary

Integrate SADIE's germline database module with existing AIRR annotation and renumbering systems. Enable germline provider selection (IMGT, OGRDB, VDJbase, custom) with backwards-compatible feature flags.

**Tech Stack**: Python 3.10+, pyhmmer, Biopython, pydantic, pytest
**Integration Points**: 3 (IgBLAST paths, HMM generation, Reference system)
**Testing Strategy**: Mirror critical path tests from existing airr/renumbering suites

## Implementation Strategy

**MVP Scope**: User Story 1 (AIRR annotation with germline provider selection)
**Delivery Order**: US1 → US2 → US3 → US4 (by priority)
**Parallel Opportunities**: US1 and US2 can be developed in parallel after foundational phase

## Prerequisites

**⚠️ External Dependency Required**: BLAST+ tools (makeblastdb, blastn)

The germlines integration requires BLAST+ to build IgBLAST databases. Install via:
- **macOS**: `brew install blast` (requires Homebrew permissions)
- **Linux**: `sudo apt-get install ncbi-blast+` or `conda install -c bioconda blast`
- **Conda**: `conda install -c bioconda blast`

Once installed, run: `python -c "from sadie.germlines import update_databases; update_databases('human')"`

## Task Progress

- **Total Tasks**: 51
- **Completed**: 34 (67%)
- **In Progress**: 0
- **Blocked**: 3 (requires gapped AA sequences for LocalHMMBuilder)
- **Remaining**: 14 (27%)

---

## Phase 1: Setup & Verification

**Goal**: Verify germlines module is ready and databases are built

- [x] T001 Verify germlines module installation at src/sadie/germlines/
- [x] T002 Verify normalized germline data exists in src/sadie/germlines/normalized/human/
- [x] T003 Build IgBLAST databases using update_databases("human") in src/sadie/germlines/
- [x] T004 Verify BLAST database files exist at src/sadie/germlines/igblast/database/human/

---

## Phase 2: Foundational Components

**Goal**: Create shared adapters and integration utilities used by all user stories

**Blocking Prerequisites**: Must complete before any user story implementation

- [x] T005 [P] Create G3 adapter class in src/sadie/germlines/g3_adapter.py
- [x] T006 [P] Implement GermlineToG3Adapter.to_g3_format() method with region extraction
- [x] T007 [P] Implement GermlineToG3Adapter.to_g3_format_batch() for bulk conversion
- [x] T008 [P] Create LocalHMMBuilder class in src/sadie/germlines/renumbering_integration.py

**Test Coverage**:
- Foundational components have unit tests in src/sadie/germlines/tests/

---

## Phase 3: User Story 1 - AIRR Annotation with Germline Provider Selection (P1)

**Story Goal**: Enable researchers to select germline database providers for AIRR annotation

**Independent Test Criteria**:
- ✅ User runs AIRR annotation with provider="imgt" → annotations reference IMGT germlines
- ✅ User runs AIRR annotation with provider="custom" → custom genes included with priority
- ✅ User runs AIRR annotation without provider parameter → default priority order used
- ✅ Existing tests in tests/unit/airr/ continue to pass (backwards compatibility)

### Tasks

- [x] T009 [P] [US1] Update GermlineData.__init__() in src/sadie/airr/igblast/germline.py with feature flag check
- [x] T010 [P] [US1] Add germlines module path logic to GermlineData class
- [x] T011 [P] [US1] Add legacy G3 path fallback logic to GermlineData class
- [x] T012 [US1] Test IgBLAST integration with feature flag enabled ✅ (BLAST+ installed, databases built)
- [x] T013 [US1] Test IgBLAST integration with feature flag disabled (backwards compat)
- [x] T014 [US1] Validate AIRR annotation results match between G3 and germlines backends

---

## Phase 4: User Story 2 - Renumbering with Germline Provider Selection (P1)

**Story Goal**: Enable researchers to select germline database providers for renumbering/HMM alignment

**Independent Test Criteria**:
- ✅ User runs renumbering with germlines backend → HMM models built from germline sequences
- ✅ User switches providers → renumbering results reflect new provider
- ✅ Existing tests in tests/unit/renumbering/ continue to pass (backwards compatibility)

### Tasks

- [x] T015 [P] [US2] Implement LocalHMMBuilder.get_hmm() in src/sadie/germlines/renumbering_integration.py
- [x] T016 [P] [US2] Implement LocalHMMBuilder._build_hmm() with Stockholm generation
- [x] T017 [P] [US2] Implement LocalHMMBuilder._get_vj_alignment_pairs() querying germlines
- [x] T018 [P] [US2] Add feature flag check function to src/sadie/renumbering/aligners/hmmer.py
- [x] T019 [US2] Update HMMER.get_hmm_models() with LocalHMMBuilder integration
- [ ] T020 [US2] Test renumbering with LocalHMMBuilder (feature flag enabled)
- [ ] T021 [US2] Test renumbering with G3 fallback (feature flag disabled)
- [ ] T022 [US2] Validate renumbering results match between G3 and germlines backends

---

## Phase 5: Reference System Integration (Foundational for US1/US2)

**Goal**: Enable Reference system to query germlines module with G3 format compatibility

**Note**: This supports both US1 and US2 by providing consistent germline data access

### Tasks

- [x] T023 [P] Add use_germlines parameter to Reference.__init__() in src/sadie/reference/reference.py
- [x] T024 [P] Initialize GermlineManager and G3Adapter in Reference.__init__() when use_germlines=True
- [x] T025 Update Reference._get_gene() to query germlines module when enabled
- [x] T026 Update Reference._get_genes() to query germlines module when enabled
- [ ] T027 Test Reference system with germlines backend (use_germlines=True)
- [ ] T028 Test Reference system with G3 backend (use_germlines=False, default)
- [ ] T029 Validate Reference output format consistency between backends

---

## Phase 6: User Story 3 - Mirrored Test Suite (P2)

**Story Goal**: Create test suite validating germlines integration produces equivalent results

**Independent Test Criteria**:
- ✅ All mirrored tests in tests/unit/germlines/ pass with germlines backend
- ✅ Test results equivalent to existing G3-based tests
- ✅ Tests run in CI alongside existing test suites

### Tasks

- [x] T030 [P] [US3] Create tests/unit/germlines/ directory
- [x] T031 [P] [US3] Create test_airr_integration.py mirroring critical AIRR tests
- [x] T032 [P] [US3] Implement test_airr_annotation_with_germlines() test case
- [x] T033 [P] [US3] Implement test_provider_selection() test case
- [x] T034 [P] [US3] Implement test_offline_operation() test case for AIRR
- [x] T035 [P] [US3] Create test_renumbering_integration.py mirroring critical renumbering tests
- [x] T036 [P] [US3] Implement test_hmmer_with_local_builder() test case
- [x] T037 [P] [US3] Implement test_hmm_caching() test case
- [x] T038 [P] [US3] Implement test_offline_operation() test case for renumbering
- [ ] T039 [US3] Run full test suite and validate all tests pass

---

## Phase 7: User Story 4 - Offline Operation (P3)

**Story Goal**: Ensure AIRR and renumbering work completely offline

**Independent Test Criteria**:
- ✅ With germlines populated and network disabled → AIRR annotation succeeds
- ✅ With cached HMMs and network disabled → renumbering succeeds
- ✅ Clear error messages if germlines not populated (first-time setup)

### Tasks

- [ ] T040 [P] [US4] Add network availability check to integration tests
- [ ] T041 [US4] Test AIRR annotation with network disabled (must succeed)
- [ ] T042 [US4] Test renumbering with network disabled (must succeed)
- [ ] T043 [US4] Test first-time setup error handling (germlines not populated)
- [ ] T044 [US4] Validate error messages are clear and actionable

---

## Phase 8: Polish & Documentation

**Goal**: Finalize documentation, performance validation, and migration guide

### Tasks

- [ ] T045 [P] Update INTEGRATION_GUIDE.md with Reference system integration details
- [ ] T046 [P] Document SADIE_USE_GERMLINES_MODULE environment variable in README.md
- [ ] T047 [P] Create migration guide for existing users in docs/
- [ ] T048 [P] Add performance benchmarking tests (germlines vs G3)
- [ ] T049 Run performance comparison and document results
- [ ] T050 Update quickstart.md with any implementation changes
- [ ] T051 Review and update API documentation for all modified modules

---

## Dependency Graph

```
Phase 1 (Setup) → Phase 2 (Foundational)
                    ↓
                    ├→ Phase 3 (US1: AIRR) ←→ Phase 5 (Reference)
                    │                          ↓
                    ├→ Phase 4 (US2: Renumbering)
                    │   ↓
                    └→ Phase 6 (US3: Tests) → Phase 7 (US4: Offline)
                        ↓
                    Phase 8 (Polish)
```

**Story Dependencies**:
- US1 and US2 are INDEPENDENT after Phase 2 (can develop in parallel)
- US3 depends on US1 and US2 completion (mirrors their functionality)
- US4 depends on US3 (uses test infrastructure)
- Phase 5 (Reference) supports both US1 and US2

---

## Parallel Execution Examples

### Phase 2 (Foundational) - All Parallel
```bash
# All adapter/builder tasks can run simultaneously
parallel_tasks=(T005 T006 T007 T008)
```

### Phase 3 (US1) - Parallel Implementation
```bash
# After foundational phase, these can run in parallel
parallel_tasks=(T009 T010 T011)
# Then sequential: T012 → T013 → T014 (testing phase)
```

### Phase 4 (US2) - Parallel Implementation
```bash
# Can run in parallel with Phase 3
parallel_tasks=(T015 T016 T017 T018)
# Then sequential: T019 → T020 → T021 → T022
```

### Phase 5 (Reference) - Parallel Implementation
```bash
# Can start after Phase 2, parallel with US1/US2
parallel_tasks=(T023 T024)
# Then sequential: T025 → T026 → T027 → T028 → T029
```

### Phase 6 (US3) - Parallel Test Creation
```bash
# All test file creation can happen in parallel
parallel_tasks=(T030 T031 T032 T033 T034 T035 T036 T037 T038)
# Then sequential: T039 (run all tests)
```

### Phase 8 (Polish) - All Parallel
```bash
# All documentation tasks can run simultaneously
parallel_tasks=(T045 T046 T047 T048)
# Then sequential: T049 → T050 → T051
```

---

## Testing Strategy

**Backwards Compatibility** (Critical):
- All existing tests must pass with feature flag disabled
- Test suites: tests/unit/airr/, tests/unit/renumbering/

**Integration Validation**:
- New tests in tests/unit/germlines/ verify equivalence
- Compare G3 vs germlines backend outputs for same inputs

**Offline Operation**:
- Network-disabled tests verify local-first operation
- Error handling tests for unpopulated databases

---

## Current Implementation Status

**Completed** (as of 2026-01-20):
- ✅ Phase 1: Setup verified (databases exist, need build step)
- ✅ Phase 2: Foundational adapters created
- ✅ Phase 3: IgBLAST integration implemented
- ✅ Phase 4: HMM integration implemented
- ✅ Phase 5: Reference system integrated
- ✅ Phase 6: Test suite created (82% - files created, APIs fixed)
- ⚠️ T013: Backwards compatibility verified (G3 fallback works)
- ❌ Phase 7-8: Offline testing and polish not started

**Blocked Tasks** (Requires BLAST+ installation):
1. T003: Build IgBLAST databases using update_databases("human")
2. T012: Test IgBLAST integration with feature flag enabled
3. T020: Test renumbering with LocalHMMBuilder

**Next Priority Tasks** (After BLAST+ installed):
1. Install BLAST+ tools (makeblastdb, blastn)
2. T003-T004: Build germlines databases
3. T012-T014: Test IgBLAST integration
4. T020-T022: Test renumbering integration
5. T027-T029: Test Reference system integration
6. T039: Run full test suite

---

## Success Criteria Validation

- **SC-001**: Mirrored AIRR tests in tests/unit/germlines/ → T031-T034
- **SC-002**: Mirrored renumbering tests in tests/unit/germlines/ → T035-T038
- **SC-003**: AIRR works with any provider → T012-T014
- **SC-004**: Renumbering works with any provider → T020-T022
- **SC-005**: Results match G3 IMGT → T014, T022, T029
- **SC-006**: Offline operation → T040-T044
- **SC-007**: No breaking changes → T013, T021, T028

---

## Notes

**Feature Flag**: `SADIE_USE_GERMLINES_MODULE` environment variable
- Default: "true" (use germlines module)
- Set to "false" for G3 API fallback (backwards compatibility)

**Performance Targets**:
- Gene lookup: <200ms (10x faster than G3)
- HMM build (uncached): <15s (equivalent to G3)
- HMM load (cached): <100ms (100x faster than G3)

**Constitution Compliance**:
- ✅ Principle V: Integration Compatibility (all gates pass)
- ✅ Principle III: Local-First Operation
- ✅ Principle I: Provider-Based Architecture

