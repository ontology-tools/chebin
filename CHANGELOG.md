# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-10-09

### Added

- InChI input. Study sets can be given as InChI as well as SMILES and ChEBI IDs,
  on the website and in the `_from_smiles` package functions, and the three can
  be mixed. An InChI is parsed into a molecule and then resolved exactly like a
  SMILES; the Chebifier API, which accepts SMILES only, is sent the canonical
  SMILES it was parsed into. A bare InChIKey is not accepted.
- ChEBI IDs are accepted in the formats they are commonly written in:
  `CHEBI:17079`, `chebi:17079`, `ChEBI:17079`, `CHEBI_17079`, `CHEBI 17079`,
  `CHEBI ID: 17079`, the bare number `17079` and the full OBO IRI. Recognition
  lives in one place (`chebin.calculations.chebi_ids`), shared by the website
  and the package.
- A controls panel on the graph canvas for hiding nodes, restoring hidden nodes,
  resetting the view and hiding labels, alongside the existing Show/Hide menu.
  Hidden nodes now stay hidden when the p-value slider is moved. Also in the
  standalone HTML export.
- Graph nodes link to their ChEBI entry: the ID in a node's tooltip is a link,
  and double-clicking a node opens the same page. The tooltip stays open long
  enough to reach the link. Also in the standalone HTML export.
- A choice of colour scale in the graph legend: **relative** (darkest = the
  smallest p-value in this graph, the previous behaviour) or **absolute**
  (darkest = a fixed cutoff, p ≤ 1e-4 by default), so colours can be compared
  between graphs. Also in the standalone HTML export.

### Changed

- Structures without a ChEBI entry of their own are classified with Chebifier's
  `best_model`, replacing the withdrawn `ELECTRA (ChEBI50-3STAR)` model.
- Chebifier-predicted parent classes are filtered by the same two rules for
  study sets as for building the restricted backgrounds: only the deepest of
  several predicted parents is kept, and a class with more than 150 leaf
  descendants is not expanded (`chebin.calculations.predicted_parents`).
  Previously a generic prediction in a study set could expand to tens of
  thousands of leaves. This can change results for SMILES/InChI study sets that
  rely on predicted parents.
- InChIKeys for ChEBI leaf classes are taken from the ChEBI OWL file rather than
  computed from each class's SMILES with RDKit, which is used only as a fallback
  for classes without one. This changes which compounds match by InChIKey, both
  for study sets and when building the backgrounds. Regenerate the data folder
  with `create_all_files` to pick it up.
- The two Homo sapiens backgrounds are now named "broad" (HMDB + LOTUS) and
  "narrow" (Recon3D) throughout the website and documentation. The short names
  passed to `narrow_background_leaves_json` (`"human"`, `"endogenous_human"`)
  are unchanged.
- Renamed the module `chebin.calculations.visualitations_and_pruning` to
  `chebin.calculations.visualizations_and_pruning`, and its function
  `graph_to_cytospace_json` to `graph_to_cytoscape_json`, to fix their spelling.
  Code that imports from the old module path must be updated; imports from the
  top-level `chebin` package are unaffected.
- The graph's node tooltip no longer uses qTip2, which closed the tooltip as
  soon as the pointer left the node. The website and the standalone export now
  share one tooltip implementation, and the graph page no longer loads jQuery,
  qTip2 or `cytoscape-qtip`.
- Removed dead code and development leftovers (commented-out code, the `old/`
  debugging scripts, a notebook and a test script inside the package).

### Fixed

- Structures RDKit cannot parse no longer make the website fail. They are
  reported on the results page, and in the `invalid_structures` list of the
  `_from_smiles` functions' diagnostics, and left out of the analysis.
- Study sets with ChEBI IDs in an unrecognised format no longer run as entities
  that match nothing, which produced an empty result with no explanation. A
  study set of bare numbers (`17079 17080`) is no longer read as one ID plus a
  weight.
- Graph menu and toolbar controls that contain an icon did nothing when the icon
  itself was clicked.
- The p-value filter's drag handle shows its grip icon again (a Font Awesome 4
  icon name that no longer exists in v6).
- Type checking in CI and locally now uses the same pinned versions, so CI no
  longer fails on diagnostics a local `prek` run reports as clean.

## [1.0.1] - 2026-09-07

### Added

- GitHub Actions workflow that publishes to PyPI when a version tag is pushed.
- `DATA_WORKFLOW.md`, a stage-by-stage description of how the data files are
  generated.

### Changed

- README revised and restructured.
- Package metadata and citation metadata updated.

## [1.0.0] - 2026-08-28

### Added

- Comprehensive type hints for all core calculation functions (Python 3.12+
  syntax)
- 173 unit tests covering enrichment analysis, visualization, and edge cases
- Integration tests validating end-to-end enrichment pipeline
- Great-docs configuration for API documentation generation
- Edge case tests for label cleaning and ChEBI ID extraction
- Support for restricted background enrichment analysis (Homo sapiens broad and
  narrow, Arabidopsis thaliana)
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

- Type checking errors (pyright) for all calculation functions
- Linting issues (ruff) in test files
- Return type mismatches in weighted enrichment functions
- Session parameter handling for None values

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

[Unreleased]: https://github.com/ontology-tools/chebin/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/ontology-tools/chebin/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/ontology-tools/chebin/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/ontology-tools/chebin/releases/tag/v1.0.0
[0.1.0]: https://github.com/Adafede/chebin/releases/tag/v0.1.0
