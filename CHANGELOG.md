# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- ChEBI IDs are accepted in every format users write them in --- `CHEBI:17079`,
  `chebi:17079`, `ChEBI:17079`, `CHEBI_17079`, `CHEBI 17079`, `CHEBI ID: 17079`,
  the bare number `17079`, and full IRIs --- rather than `CHEBI:17079` alone.
  Recognition lives in one place (`chebin.calculations.chebi_ids`), shared by
  the website's study-set parser, `normalize_id()` and the SMILES/ChEBI-ID check
- Comprehensive type hints for all core calculation functions (Python 3.12+
  syntax)
- 173+ unit tests covering enrichment analysis, visualization, and edge cases
- Integration tests validating end-to-end enrichment pipeline
- Great-docs configuration for API documentation generation
- Edge case tests for label cleaning and ChEBI ID extraction
- Support for narrow background enrichment analysis (human, A. thaliana, ReconX)
- Multiple testing correction: Bonferroni and Benjamini-Hochberg FDR methods
- Graph pruning strategies: root-children, linear-branch, high-p-value,
  zero-degree
- Visualization utilities: graph construction, node filtering, ID normalization

### Changed

- Fixed return type annotations for functions returning 4-tuples
- Improved website enrichment endpoint to handle classification parameter
  properly
- Enhanced error handling in edge cases
- Refactored test structure for better organization

### Fixed

- Study sets written with an unexpected ChEBI ID format no longer run through as
  entities that match nothing, which produced an empty result with no
  explanation. A bare-number study set (`17079, 17080`) is also no longer read
  as one ID plus a weight
- Type checking errors (pyright) for all calculation functions
- Linting issues (ruff) in test files
- Return type mismatches in weighted enrichment functions
- Session parameter handling for None values

### Testing

- Added 27 visualization and pruning strategy tests
- Added 24 integration tests for enrichment pipeline
- Added 9 edge case tests for label utilities
- Total: 173 passing tests across 8 test modules

## [0.1.0] - 2024-08-10

### Added

- Initial project structure
- ChEBI ontology enrichment analysis tool
- Fisher's exact test implementation for enrichment analysis
- Weighted enrichment analysis using Lugannani-Rice saddlepoint approximation
- Graph-based visualization and pruning of enrichment results
- Web interface for enrichment analysis
- Narrow background support for species-specific analysis
- ChEBI ID and label utilities
- Pre-Fisher's calculations for data preparation

### Features

#### Core Calculations

- `calculate_p_value()`: Fisher's exact test for 2x2 contingency tables
- `run_enrichment_analysis()`: Full enrichment pipeline with multiple strategies
- `run_enrichment_analysis_plain_enrich_pruning_strategy()`: Fixed pruning
  strategy
- `calculate_weighted_p_value()`: Weighted variant of Fisher's test
- `run_weighted_enrichment_analysis()`: Full weighted enrichment pipeline
- Multiple testing correction (Bonferroni, Benjamini-Hochberg FDR)

#### Visualization & Pruning

- Graph construction from enrichment results
- Root-children pruning (distance from root)
- Linear-branch pruning (collapse short branches)
- High-p-value pruning (remove weak enrichments)
- Zero-degree pruning (remove isolated nodes)
- Graph-to-JSON export for visualization

#### Data Processing

- ChEBI data loading and normalization
- ID-to-name mapping
- Role classification (structural, functional)
- Leaf and ancestor computation
- Narrow background leaf restriction

### Dependencies

- Python 3.12+
- networkx: Graph algorithms
- rdkit: Chemical structure processing
- requests: HTTP client
- flask: Web framework
- pyhornedowl: OWL ontology parsing

### Quality Assurance

- Pre-commit hooks via prek
- Type checking with pyright
- Code linting with ruff
- Formatting with black and ruff-format
- Dependency auditing with deptry
- Test framework: pytest

--------------------------------------------------------------------------------

[Unreleased]: https://github.com/Adafede/chebin/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Adafede/chebin/releases/tag/v0.1.0
