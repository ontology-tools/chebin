# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- The graph page has a controls panel sitting on the canvas itself, so hiding
  nodes, restoring the original view and hiding labels no longer have to be
  found in the navbar's Show/Hide menu. It mirrors the p-value filter on the
  opposite corner --- same chrome, same drag handle --- and reports how many
  nodes are hidden. The menu keeps its items; both routes run the same code and
  stay in step, so toggling labels from one updates the other. "Hide selected
  nodes" also works on a single node from its right-click menu, and because the
  selection tools already offer "select descendants" a whole subtree goes in two
  clicks. Hidden nodes now stay hidden: the action used to be undone by the next
  nudge of the threshold slider, which made it look as though nothing had
  happened. Two ways to bring them back: "Return hidden nodes" restores what is
  hidden and nothing else, leaving the layout, labels, zoom and threshold as
  they are, so nodes reappear only if they still pass the p-value filter; "Reset
  graph" goes further and returns the whole view to how it opened, layout and
  zoom included. Nothing is deleted from the underlying graph, so the colour
  scale does not shift as nodes go. The standalone export gets the panel too

- Graph nodes link to their ChEBI entry. A node's hover tooltip is titled
  `D-glucoside (CHEBI:35436)` with the id itself a link to that entry, and
  double-clicking a node opens the same page directly. Links go to
  `https://www.ebi.ac.uk/chebi/CHEBI:35436`, ChEBI's current entry page --- the
  older `searchId.do` and `chebiOntology.do` forms now only redirect there. The
  tooltip is parked flush against the node rather than following the cursor, so
  the link is something to aim at rather than chase, and it survives the pointer
  leaving the node by 600 ms --- indefinitely once the pointer is on the tooltip
  itself. Both the website and the offline standalone export have it; in the
  export the link is inert until clicked, so the page still needs no network to
  render. Nodes carry their id as a `chebi_id` CURIE for this, and graphs
  written before that field fall back to the node's OBO IRI, so already exported
  pages link too

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

- The graph's node tooltip no longer uses qTip2. The `cytoscape-qtip` extension
  binds `hide.event` as a Cytoscape event and calls `qtipApi.hide()` itself, so
  qTip2's own `hide.fixed` and `hide.delay` never ran: the tooltip closed the
  instant the pointer left the node, whatever delay was configured, and any pan
  or zoom dismissed it too. That left the ChEBI link in it unreachable. It is
  now the same hand-rolled tooltip the standalone export already used, so both
  views share one implementation instead of two that drift. This drops the
  jQuery, qTip2 and `cytoscape-qtip` CDN dependencies from the graph page
- `calculate_weighted_p_value` annotates its `saddler` parameter as `_SaddleSum`
  instead of `object`, so the `saddler.pvalue()` call is type checked rather
  than silenced by a blanket `# type: ignore`
- Dropped a redundant truthiness check on `ET.ParseError.position` in the HMDB
  extractor. A non-empty tuple is always truthy, so the guard never varied
- Fixed return type annotations for functions returning 4-tuples
- Improved website enrichment endpoint to handle classification parameter
  properly
- Enhanced error handling in edge cases
- Refactored test structure for better organization

### Fixed

- The p-value filter's drag handle shows its grip icon again. It asked for
  `fa fa-arrows`, a Font Awesome 4 name that was renamed in v6 and has no
  compatibility shim loaded, so the icon had been rendering as nothing

- Graph menu and toolbar controls read their action from the element the handler
  is bound to rather than from whatever was clicked. Every control was
  text-only, so the two were the same thing until a control contained an icon,
  at which point clicking the icon read the icon's own (absent) action and did
  nothing

- Study sets written with an unexpected ChEBI ID format no longer run through as
  entities that match nothing, which produced an empty result with no
  explanation. A bare-number study set (`17079, 17080`) is also no longer read
  as one ID plus a weight

- CI type-check failures that could not be reproduced locally. The `ty-check`
  hook in `prek.toml` left its dependencies unpinned, so CI resolved the latest
  `ty` on every run while the local hook environment stayed cached at an older
  one --- CI failed on diagnostics a local `prek` run reported as clean. `ty`
  and the libraries it reads as type information are now pinned to the versions
  in `uv.lock`, so local and CI resolve identically

- `rdkit` added to the `ty-check` hook's dependencies. The hook resolves imports
  from its own environment rather than the project `.venv`, so without it every
  `rdkit` import was an unresolved-import masked by a blanket `# type: ignore`

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
