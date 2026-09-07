# Datafiles Workflow

This file explains how the generated data files are built, stage by stage. It is
reference material for contributors regenerating the data folder or changing the
pipeline --- you do **not** need any of it to use the `chebin` package or the
web application.

See the [main README](README.md) for installation and package usage.

In the web application, the files are regenerated automatically once a month.

<details>
<summary><strong>0. Python environment</strong> --- only needed if you are
working from the repository.</summary>

Dependencies and the Python version requirement (`>=3.12`) are declared in
`pyproject.toml`, and `uv.lock` pins every dependency (and transitive
dependency) to an exact version, so the environment is fully reproducible across
machines.

To create it, install [uv](https://docs.astral.sh/uv/) and run, from the
repository root:

```bash
uv sync
```

This creates (or updates) `.venv/` with every package pinned to the exact
version in `uv.lock`. `.python-version` pins the interpreter itself (3.14);
`uv sync` downloads a matching Python automatically if one isn't already
available, so no separate Python install step is needed.

Development-only tools (`pytest`, `prek` linting, `coverage`, `great-docs`) are
declared as a separate dependency group and are included by the plain `uv sync`
above. To skip them, e.g. for a production-only install, use `uv sync --no-dev`.

If you'd rather not use uv, the runtime dependencies (see `pyproject.toml` for
the authoritative, version-constrained list) can be installed manually with:

`pip install flask networkx numpy pandas py-horned-owl rdkit requests scipy`

</details>

Note that the steps below do not need to be run individually: the files can all
be built through `jobs/run_create_files.sh` or the package's
`create_all_files()`. This breakdown is here for understanding what each stage
does.

## 1. Load ChEBI

The ChEBI ontology is downloaded from
https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl by
`src/chebin/preparing_data/load_chebi.py` (task *"download_and_load_chebi"*) and
saved as `data/chebi.owl`. It is re-downloaded on every run by default, since
ChEBI is updated continuously; pass `force=False` to reuse an existing local
copy instead.

## 2. Remove leaf classes and save maps

**What counts as a leaf:** A class is a leaf if it has its own valid SMILES
string that is *not* a wildcard/R-group placeholder (SMILES containing the dummy
atom `*`, e.g. `*C(N)C(=O)O`, are rejected via RDKit). This holds **regardless
of whether the class has subclasses** --- a parent class with a proper SMILES
(e.g. `proline`) is a leaf in its own right, alongside its more specific
children (e.g. `L-proline`, `D-proline`).

**Flattening the hierarchy (the "splice"):** Because a leaf must be terminal,
whenever a class sits under a leaf it is reconnected to that leaf's nearest
*non-leaf* ancestors, climbing past chains of stacked leaves. So `L-proline`
stops pointing at the leaf `proline` and instead points directly at the real
category above it (e.g. `alpha-amino acid`), while keeping every other ancestor
it already had (e.g. `D-proline` keeps `D-alpha-amino acid`). After the splice,
no leaf is any class's parent, so all SMILES-bearing classes become siblings
under the genuine (non-leaf) category terms.

Run task *"remove_leaves_with_smiles"* to find leaf classes, splice the
hierarchy, and filter out deprecated classes. The following files are created:

- A filtered OWL file with the remaining classes (from `save_filtered_owl`). The
  current file in this workspace is
  `data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl`.
- A **flattened** subclass map JSON file (`data/chebi_subclass_map.json`)
  mapping all classes to their direct subclasses after the splice (from
  `find_leaf_classes_with_smiles_and_deprecated`). Leaves have no subclasses
  here; the file also includes deprecated classes.
- A leaf-to-parents map JSON file mapping each leaf class to all of its
  **non-leaf** ancestors (from `find_leaf_classes_with_smiles_and_deprecated`).
  This file is used in later calculations.

Run task *"build_parent_map"* to create:

- `data/chebi_parent_map.json`, a map of all classes to their direct parents
  after the splice (deprecated classes are excluded). It is derived from the
  flattened subclass map, so the graph built from it shows the same sibling
  structure.

Run task *"map_names_to_classes"* to build:

- `data/chebi_id_to_name_map.json`, which maps short CHEBI IDs (e.g.
  `CHEBI_111`) to their names.

## 2.5 Save maps connected to the roles of the classes

Maps that include the roles of the classes are needed for some enrichment
calculations. These are made in
`src/chebin/calculations/prepare_role_calculations.py`.

First, run the task *"find has_role connections"*. This parses the OWL file
directly and produces a map from all classes to their **direct** roles (not
including any roles that ancestors have):

- `data/class_to_direct_roles_map.json`

This task also calls `create_leaves_to_all_roles_map`, which writes

- `data/removed_leaf_classes_to_ALL_roles_map.json` (using
  `data/removed_leaf_classes_to_ALL_parents_map.json` and
  `data/chebi_parent_map.json`)

Second, run the task *"build leaf to all roles map"*. This builds:

- `data/removed_leaf_classes_to_ALL_roles_map.json`
- `data/roles_to_leaves_map.json`

The first file maps each removed leaf class to (a) its direct roles, (b) roles
inherited from ancestor classes, and (c) **ancestors of those roles** in the
role hierarchy. The second file is the inverse: it maps each role class to all
leaf classes connected to that role.

Note: because the hierarchy was flattened in step 2, a leaf inherits roles only
from its **non-leaf** ancestors. It no longer inherits roles asserted directly
on a leaf-ancestor --- e.g. `D-proline` no longer inherits roles (such as *human
metabolite*) that are asserted on the generic `proline`, which is now itself a
leaf. A leaf's own direct roles are always kept.

Third, run the task *"build class to all roles map"* to create:

- `data/class_to_all_roles_map.json`

Here, each class is mapped to its direct roles, roles inherited from ancestor
classes, and **descendants** of those roles in the role hierarchy.

## 3. Split up the ontology based on structure

Classes are sorted into structural vs. functional (role) sets by
`identify_structural_vs_functional()` in
`src/chebin/preparing_data/pruning_split_up_structure.py`, which walks the
descendants of the structural and role root classes and returns three sets of
class IRIs: `structural_classes`, `functional_classes`, and `unknown_classes`
(classes under neither root, kept only for troubleshooting). These three sets
feed the `Classification` column of the CSV created in step 4 below.

There are two ways to obtain these sets, and
`src/chebin/preparing_data/create_files.py` uses the fast one:

- **Automated (fast) path --- used by
  `src/chebin/preparing_data/create_files.py`:**
  `identify_structural_vs_functional()` is called with the already-in-memory
  flattened subclass map (`data/chebi_subclass_map.json`, built in step 2), so
  descendants are found via plain dict lookups instead of per-node ontology API
  calls. The resulting class sets are passed straight into
  `save_leaf_classes_with_smiles()` (step 4) without ever touching disk --- no
  OWL files are written for this step.
- **Manual path --- `src/chebin/preparing_data/pruning_split_up_structure.py`,
  task *"split_structural_functional"*:** Run standalone, the function is called
  without the subclass map, so it falls back to the slower ontology API. The
  three class sets are then written out as separate OWL files via
  `split_owl_by_type()`, creating `_structural.owl`, `_functional.owl`, and
  `_unknown.owl` versions of the previously filtered ontology (e.g.
  `data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl`). The
  *"save_removed_leaf_classes"* task (step 4) in
  `src/chebin/preparing_data/pruning_smiles.py` can then load these OWL files
  back in to recover the same three class sets, if run outside of
  `src/chebin/preparing_data/create_files.py`.

Since the splice in step 2 preserves every class's reachable non-leaf ancestors,
descendant sets for these (non-leaf) roots are the same whether collected from
the flattened subclass map or from live ontology traversal --- so the two paths
produce identical class sets. The OWL files are just an on-disk,
human-inspectable form of the same information, useful for manual runs or
debugging.

## 4. Save a file with the removed leaf classes

Go back to `src/chebin/preparing_data/pruning_smiles.py` and run task
*"save_removed_leaf_classes"* to save the removed leaf classes in a CSV file.
The current file in this workspace is

- `data/removed_leaf_classes_with_smiles.csv`

The CSV contains `IRI`, `SMILES`, and `Classification`, where the classification
is inferred from the class's direct parents in the structural/functional split.
Every leaf has a row here, including classes that have subclasses but carry
their own valid SMILES (e.g. `proline`); classes whose SMILES is a
wildcard/R-group placeholder are not leaves and do not appear. In ChEBI, classes
with SMILES are expected to fall under structural roots, so entries classified
as **functional** are likely misclassified and are excluded from downstream
calculations (they are kept in the file for troubleshooting).

## 5. Fisher's Calculations

First (only needed once), run task *"build_class_to_leaf_map"* in
`src/chebin/calculations/pre_fishers_calculations.py` to create
`data/class_to_leaf_descendants_map.json`, which maps each class to all of its
removed leaf descendants using
`data/removed_leaf_classes_to_ALL_parents_map.json`. Leaf classes are **not**
keys in this map: after the splice no leaf appears as another class's ancestor,
so only the genuine (non-leaf) category terms become keys. A leaf is therefore
only ever counted as a member of its categories, never tested as a category
itself.

Enrichment calculations can be run in
`src/chebin/calculations/fishers_calculations.py`, but this is most easily done
via the web application. Either use the website link (easiest since no
preparation steps to obtain all the necessary files are needed) or run
`website/app.py` locally.

## 6. Needed for human dataset

1. Download LOTUS compound--taxon data from Wikidata via the QLever SPARQL
   endpoint using `src/chebin/preparing_data/wikidata/get_lotus.py`. This is run
   automatically by `src/chebin/preparing_data/create_files.py`, but can also be
   run standalone:

   ```bash
   python -m chebin.preparing_data.wikidata.get_lotus
   ```

   Output:

   - `data/lotus_homo_sapiens.csv`
   - `data/lotus_arabidopsis_thaliana.csv`

2. Connect the LOTUS CSVs to ChEBI IDs using `connect_lotus_csv_to_chebi_ids()`
   in `src/chebin/preparing_data/wikidata/get_wikidata_lotus.py`.

   Output:

   - `data/wikidata/created/lotus_homo_sapiens_with_chebi_ids.tsv`
   - `data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids.tsv`

3. Extract HMDB compounds using `extract_hmdb_to_file()` in
   `src/chebin/preparing_data/hmdb/extract_hmdb.py`.

   `data/hmdb_metabolites.xml` is required and must be downloaded manually from
   https://hmdb.ca/downloads (use the 'All Metabolites' XML).

   Output: `data/hmdb_metabolites_extract.tsv`

4. Filter HMDB to only keep compounds with status "quantified" or "detected"
   using `filter_hmdb_statuses_main()` in
   `src/chebin/preparing_data/hmdb/filter_hmdb_statuses.py`.

   Output: `data/hmdb_metabolites_extract_quantified_detected.tsv`

5. Find missing ChEBI IDs using `run_find_missing_chebis(source)` in
   `src/chebin/preparing_data/wikidata/find_missing_chebis.py` (also runnable
   via `jobs/run_find_missing_chebis.sh [source]`). The `source` argument must
   be one of the presets in `SOURCE_PRESETS`: `"lotus_hs"`, `"lotus_at"`, or
   `"hmdb"`.

   ChEBI ID matching is attempted in this order:

   - direct ChEBI matches (LOTUS)
   - exact SMILES match against ChEBI leaf classes
   - InChIKey match against ChEBI leaf classes
   - the Chebifier API

   Output (depending on source):

   - `data/wikidata/created/lotus_homo_sapiens_with_chebi_ids_updatedchebis.tsv`
   - `data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv`
   - `data/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv`

6. Combine HMDB and LOTUS Homo sapiens sources using `combine_datasets()` in
   `src/chebin/preparing_data/wikidata/combine_human_datasets.py`. Rows with no
   ChEBI ID are dropped.

   Output: `data/combined_hmdb_wikidata.tsv`

7. Create a file with the human leaf classes using `gather_narrow_leaves()` in
   `src/chebin/preparing_data/wikidata/narrow_background_fishers.py`.

   Files needed:

   - `compounds_tsv = "data/combined_hmdb_wikidata.tsv"`
   - `leaves_csv = "data/removed_leaf_classes_with_smiles.csv"`
   - `class_to_leaf_map = "data/class_to_leaf_descendants_map.json"`
   - `taxon_label = "homo_sapiens"` (recorded in the output JSON for
     traceability)

   Output: `data/human_entities_leaves.json`

## Narrow background for a single Wikidata taxon (e.g. Arabidopsis thaliana)

This is the same workflow as above, but since there is only one source
(Wikidata), steps 3, 4, and 6 (HMDB extraction/filtering and combining datasets)
are skipped entirely.

1. The Arabidopsis thaliana LOTUS CSV (`data/lotus_arabidopsis_thaliana.csv`)
   and its ChEBI-matched TSV
   (`data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids.tsv`) are
   already produced in steps 1--2 above.

2. Fill in any still-missing ChEBI IDs using the `"lotus_at"` preset in
   `src/chebin/preparing_data/wikidata/find_missing_chebis.py`.

   Output:
   `data/wikidata/created/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv`

3. Build the leaf classes with `gather_narrow_leaves()` in
   `src/chebin/preparing_data/wikidata/narrow_background_fishers.py`, passing
   the file from step 2 as `compounds_tsv` and
   `taxon_label="arabidopsis_thaliana"`.

Output: `data/arabidopsis_thaliana_leaves.json`

## Endogenous human background (Recon3D)

This path is independent of step 6 and the Wikidata/HMDB workflow above; it only
needs the files from steps 1--5 (`data/removed_leaf_classes_with_smiles.csv` and
`data/class_to_leaf_descendants_map.json`).

Run `src/chebin/preparing_data/BiGG/get_model.py`. This:

1. Downloads the Recon3D model JSON from BiGG
   (`http://bigg.ucsd.edu/static/models/Recon3D.json`) to `data/Recon3D.json`.
2. Calls `gather_recon3d_leaves()`, which collapses compartment-specific
   metabolite entries into unique compounds and resolves each to leaf ChEBI
   classes (directly, via parent expansion, or via UniChem InChIKey/HMDB
   cross-reference, as described above). Running the script prints a breakdown
   of how many compounds were resolved by each method, and how many were left
   unresolved.

Output: `data/recon3d_leaves.json` (same `narrow_leaves` JSON shape as the other
narrow backgrounds above, so it plugs into the website as
`NARROW_BACKGROUND_LEAVES_JSON['endogenous_human']` without further changes).
