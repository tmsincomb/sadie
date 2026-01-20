# Feature Specification: Germline Database Integration

**Feature Branch**: `002-germline-integration`  
**Created**: 2026-01-19  
**Status**: Draft  
**Input**: Connect SADIE's new germline database to existing sadie.airr.Airr and sadie.numbering.Renumbering to allow picking new germlines database, plus create tests in tests/unit/germlines mirroring airr and renumbering tests using backend germline instead of default IMGT from G3

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Select Germline Provider for AIRR Analysis (Priority: P1)

As a researcher running immune repertoire analysis, I want to select which germline database (IMGT, OGRDB, VDJbase, or custom) my AIRR annotation uses, so that I can leverage alternative gene databases or my own custom germlines instead of being limited to the default IMGT from G3.

**Why this priority**: This is the core integration point - enabling germline provider selection in the primary AIRR analysis workflow unlocks all downstream functionality.

**Independent Test**: Can be tested by running AIRR annotation with different provider parameters and verifying gene calls come from the specified database source.

**Acceptance Scenarios**:

1. **Given** a user has sequence data and the new germlines module is populated, **When** they run AIRR analysis with `provider="imgt"`, **Then** gene annotations reference IMGT-sourced germlines from the local database.

2. **Given** a user has custom germline sequences in `sources/custom/`, **When** they run AIRR analysis with `provider="custom"`, **Then** gene annotations include their custom genes with priority over default providers.

3. **Given** a user does not specify a provider, **When** they run AIRR analysis, **Then** the system uses the default priority order (custom > ogrdb > vdjbase > imgt) from the germlines module.

---

### User Story 2 - Select Germline Provider for Renumbering (Priority: P1)

As a researcher performing antibody numbering, I want to select which germline database my renumbering/HMM alignment uses, so that my numbering schemes align against the same germlines used in my analysis pipeline.

**Why this priority**: Renumbering is a parallel critical path alongside AIRR - both need germline integration for a complete solution.

**Independent Test**: Can be tested by running renumbering with different provider parameters and verifying HMM alignments use the specified germline source.

**Acceptance Scenarios**:

1. **Given** a user has antibody sequences, **When** they run renumbering with the germlines backend, **Then** HMM models are built from the selected germline provider's sequences.

2. **Given** a user switches germline providers between analyses, **When** they run renumbering, **Then** results reflect the currently selected provider's gene database.

---

### User Story 3 - Consistent Test Suite Using Germlines Backend (Priority: P2)

As a developer maintaining SADIE, I want a test suite in `tests/unit/germlines/` that mirrors the existing airr and renumbering tests but uses the local germlines backend, so that I can validate the germlines integration produces equivalent results to the G3-based tests.

**Why this priority**: Test coverage ensures the integration is correct and provides regression protection for future changes.

**Independent Test**: Can be tested by running the new test suite and comparing results with existing G3-based tests.

**Acceptance Scenarios**:

1. **Given** the existing `tests/unit/airr/` test suite, **When** I run the mirrored tests in `tests/unit/germlines/`, **Then** all tests pass using the local germlines backend instead of G3.

2. **Given** the existing `tests/unit/renumbering/` test suite, **When** I run the mirrored renumbering tests with germlines backend, **Then** all tests pass with equivalent results.

---

### User Story 4 - Offline Germline Operation (Priority: P3)

As a researcher working in an environment without reliable internet, I want AIRR and renumbering to work completely offline using the local germlines database, so that I am not dependent on G3 API availability.

**Why this priority**: Offline capability is a key benefit of the germlines module but depends on P1/P2 being complete first.

**Independent Test**: Can be tested by disabling network access and verifying all analysis workflows complete successfully.

**Acceptance Scenarios**:

1. **Given** the germlines module is populated and network is unavailable, **When** I run AIRR analysis, **Then** annotation completes successfully using local data.

2. **Given** cached HMM models exist, **When** I run renumbering offline, **Then** numbering completes without network errors.

---

### Edge Cases

- What happens when a requested germline provider has no data for the specified species?
- How does the system handle conflicting gene names across multiple providers?
- What happens when switching providers mid-analysis for linked datasets?
- How does the system behave when custom germlines have invalid sequences?
- What happens when the germlines module is not yet populated (first-time setup)?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST allow users to specify a germline provider (e.g., "imgt", "ogrdb", "vdjbase", "custom") when initializing AIRR analysis.

- **FR-002**: System MUST allow users to specify a germline provider when initializing Renumbering operations.

- **FR-003**: System MUST use a default provider priority order (custom > ogrdb > vdjbase > imgt) when no provider is explicitly specified.

- **FR-004**: System MUST use local germlines module data instead of G3 API when germlines backend is selected.

- **FR-005**: System MUST provide clear error messages when a requested provider has no data for the specified species.

- **FR-006**: System MUST support backwards compatibility - existing code using default settings continues to work without changes.

- **FR-007**: System MUST include a new test directory `tests/unit/germlines/` containing tests that mirror `tests/unit/airr/` functionality using germlines backend.

- **FR-008**: System MUST include tests that mirror `tests/unit/renumbering/` functionality using germlines backend.

- **FR-009**: System MUST support the same species, chains, and segments currently supported by the existing AIRR and Renumbering modules.

- **FR-010**: System MUST maintain identical output format/schema regardless of which germline provider is used.

### Key Entities

- **GermlineProvider**: Represents a source of germline data (IMGT, OGRDB, VDJbase, custom) with its priority level and available species.

- **GermlineSelection**: Configuration specifying which provider(s) to use for an analysis run, including priority ordering.

- **GermlineData**: The germline gene information used by IgBLAST including V, D, J segment databases and auxiliary files - now sourced from local germlines module.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of existing `tests/unit/airr/` test functionality is mirrored in `tests/unit/germlines/` with passing results using germlines backend.

- **SC-002**: 100% of existing `tests/unit/renumbering/` test functionality is mirrored with passing results using germlines backend.

- **SC-003**: Users can complete AIRR analysis with any available germline provider without errors.

- **SC-004**: Users can complete renumbering with any available germline provider without errors.

- **SC-005**: Gene annotation results using IMGT provider from germlines module match results from G3-based IMGT for the same input sequences.

- **SC-006**: System operates fully offline after initial germlines database population.

- **SC-007**: No breaking changes to existing user workflows - all existing tests in `tests/unit/airr/` and `tests/unit/renumbering/` continue to pass.

## Assumptions

- The germlines module (`src/sadie/germlines/`) is already implemented and populated with IMGT data as per the `001-germline-completion` feature.
- The INTEGRATION_GUIDE.md in the germlines module provides accurate integration patterns to follow.
- Existing test fixtures are compatible with the germlines backend or can be adapted.
- The priority ordering (custom > ogrdb > vdjbase > imgt) is the correct default based on project requirements.
