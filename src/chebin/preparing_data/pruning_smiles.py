"""Filter the ontology to remove leaf classes with SMILES annotations."""

import csv
import json
import os
import random
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict

import pandas as pd

from chebin.preparing_data.load_chebi import load_chebi, load_ontology


def _has_wildcard(smiles):
    """True if the SMILES contains a wildcard/dummy atom (ChEBI's R-group placeholder '*').

    In SMILES syntax '*' is only ever the dummy atom, so this substring test is exactly
    equivalent to checking for an atomic-number-0 atom (verified: zero false positives across
    all ChEBI SMILES) without parsing. Unlike a full RDKit parse, this keeps real compounds
    whose curated SMILES merely fail strict sanitization (kekulization/valence quirks).
    """
    return "*" in smiles


ancestor_cache = {}


def get_all_ancestors(chebi_ontology, cls, visited=None):
    """Recursively collect all ancestor classes, using cache to avoid recomputation."""
    cls_str = str(cls)
    if cls_str in ancestor_cache:
        return ancestor_cache[cls_str]

    if visited is None:
        visited = set()

    ancestors = []
    direct_parents = chebi_ontology.get_superclasses(cls)

    for parent in direct_parents:
        parent_str = str(parent)
        if parent_str not in visited:  # Avoid cycles
            visited.add(parent_str)
            ancestors.append(parent_str)
            # Recursively get ancestors of this parent
            ancestors.extend(get_all_ancestors(chebi_ontology, parent, visited))

    ancestors = list(set(ancestors))  # Remove duplicates
    ancestor_cache[cls_str] = ancestors  # Cache the result

    return ancestors


def _invert_subclass_map(subclass_map, all_class_ids=None):
    """Convert parent->children map into child->parents map."""
    parent_map = defaultdict(list)

    for parent, children in subclass_map.items():
        for child in children:
            parent_map[child].append(parent)

    # Keep all known classes present, even if they have no parents.
    if all_class_ids:
        for cls in all_class_ids:
            parent_map.setdefault(cls, [])

    # Deduplicate parent lists
    for child in parent_map:
        parent_map[child] = list(set(parent_map[child]))

    return dict(parent_map)


def _get_all_ancestors_from_parent_map(class_id, parent_map, cache):
    """Compute all ancestors using an iterative traversal over parent_map."""
    if class_id in cache:
        return cache[class_id]

    ancestors = set()
    stack = list(parent_map.get(class_id, []))

    while stack:
        parent = stack.pop()
        if parent in ancestors:
            continue
        ancestors.add(parent)

        if parent in cache:
            ancestors.update(cache[parent])
            continue

        stack.extend(parent_map.get(parent, []))

    ordered = list(ancestors)
    cache[class_id] = ordered
    return ordered


def _nearest_nonleaf_parents_map(raw_parent_map, leaf_set):
    """For every class, compute its nearest non-leaf ancestors, climbing past any leaf
    parents (handles chains of stacked leaves). Returns a child -> parents dict in which
    no leaf ever appears as a parent value."""
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))
    memo = {}

    def compute(node, stack):
        if node in memo:
            return memo[node]
        if node in stack:  # cycle guard (ChEBI is a DAG, but stay safe)
            return set()
        result = set()
        for p in raw_parent_map.get(node, ()):
            if p not in leaf_set:
                result.add(p)  # real (non-leaf) category: stop climbing
            else:
                result |= compute(p, stack | {node})  # leaf parent: climb past it
        memo[node] = result
        return result

    return {node: list(compute(node, frozenset())) for node in raw_parent_map}


def _splice_leaves_from_hierarchy(raw_parent_map, leaf_set):
    """Flatten the hierarchy so leaves are terminal: each class reconnects to its nearest
    non-leaf ancestors and leaves end up with no children. Returns
    (flattened_subclass_map, flattened_parent_map)."""
    flattened_parent_map = _nearest_nonleaf_parents_map(raw_parent_map, leaf_set)

    flattened_subclass_map = defaultdict(list)
    for child, parents in flattened_parent_map.items():
        for parent in parents:
            flattened_subclass_map[parent].append(child)
    flattened_subclass_map = {
        p: sorted(set(kids)) for p, kids in flattened_subclass_map.items()
    }
    return flattened_subclass_map, flattened_parent_map


def _verify_splice(raw_parent_map, flattened_parent_map, leaf_set, sample_size=8000):
    """Safety net run before any derived file is written. Asserts the splice preserved
    every node's reachable NON-leaf ancestors (no orphaning, no invented connections) and
    that no leaf is left as anyone's parent."""

    def closure(node, pm):
        seen = set()
        stack = list(pm.get(node, ()))
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(pm.get(x, ()))
        return seen

    nodes = list(raw_parent_map)
    random.seed(0)
    sample = nodes if len(nodes) <= sample_size else random.sample(nodes, sample_size)
    for n in sample:
        raw_nonleaf = {a for a in closure(n, raw_parent_map) if a not in leaf_set}
        spliced = closure(n, flattened_parent_map)
        assert raw_nonleaf == spliced, (
            f"Splice changed reachable non-leaf ancestors for {n}"
        )

    for child, parents in flattened_parent_map.items():
        assert not any(p in leaf_set for p in parents), (
            f"Splice left a leaf as a parent of {child}"
        )

    print(
        f"Splice verification passed: {len(sample)} sampled nodes, non-leaf ancestor sets preserved.",
    )


def _build_subclass_map_from_axioms(chebi_ontology, all_classes):
    """Build parent -> direct subclasses map in one pass over ontology axioms."""
    class_set = {str(cls) for cls in all_classes}
    subclass_map = defaultdict(list)

    for idx, axiom in enumerate(chebi_ontology.get_axioms()):
        component = axiom.component
        if type(component).__name__ != "SubClassOf":
            continue

        sub = str(component.sub)
        sup = str(component.sup)

        # Keep only named-class relations (skip restrictions/anonymous expressions)
        if sub in class_set and sup in class_set:
            subclass_map[sup].append(sub)

        if (idx + 1) % 500000 == 0:
            print(f"Processed {idx + 1} axioms for subclass mapping...")

    for parent in subclass_map:
        subclass_map[parent] = list(set(subclass_map[parent]))

    return dict(subclass_map)


def find_leaf_classes_with_smiles_and_deprecated(
    chebi_ontology,
    smiles_property,
    deprecated_property,
    subclass_map_file,
    leaf_parents_map_file,
    use_found_leaf_classes=False,
    removed_leaf_classes_file=None,
):

    print("\nFinding leaf classes with SMILES...")

    # Get all classes
    all_classes = chebi_ontology.get_classes()
    all_class_ids = [str(cls) for cls in all_classes]
    print(f"Total classes: {len(all_classes)}")

    # Load or build subclass mapping for faster lookup
    try:
        subclass_map_start = time.time()
        with open(subclass_map_file) as f:
            subclass_map = json.load(f)
        subclass_map_elapsed = time.time() - subclass_map_start
        print(
            f"Loaded subclass mapping from {subclass_map_file} with {len(subclass_map)} entries "
            f"in {subclass_map_elapsed:.2f} seconds.",
        )
    except FileNotFoundError:
        print(
            "Building subclass mapping from axioms in a single pass (may take a while)...",
        )
        subclass_map_start = time.time()

        compute_start = time.time()
        subclass_map = _build_subclass_map_from_axioms(chebi_ontology, all_classes)
        compute_elapsed = time.time() - compute_start

        write_start = time.time()
        with open(subclass_map_file, "w") as f:
            json.dump(subclass_map, f)
        write_elapsed = time.time() - write_start

        subclass_map_elapsed = time.time() - subclass_map_start
        print(
            f"Saved subclass mapping to {subclass_map_file} in {subclass_map_elapsed:.2f} seconds "
            f"({subclass_map_elapsed / 60:.2f} minutes).",
        )
        print(
            f"Subclass map timing breakdown: compute={compute_elapsed:.2f}s, "
            f"json_write={write_elapsed:.2f}s.",
        )

    # Build child -> direct parents map once and reuse it for fast ancestor lookup.
    direct_parent_map = _invert_subclass_map(subclass_map, all_class_ids)

    # Precompute property strings
    smiles_prop_str = f"<{smiles_property}>"
    deprecated_prop_str = f"<{deprecated_property}>"

    classes_with_smiles = []
    leaf_classes_with_smiles = []
    deprecated_classes = []

    if use_found_leaf_classes and removed_leaf_classes_file:
        print(
            f"Loading previously found leaf classes with SMILES from {removed_leaf_classes_file}...",
        )
        df = pd.read_csv(removed_leaf_classes_file)
        leaf_classes_with_smiles = df["IRI"].tolist()
        print(
            f"Loaded {len(leaf_classes_with_smiles)} leaf classes with SMILES from file.",
        )

    else:  # Scan ontology to find leaf classes with SMILES
        # Counters
        i = 0
        j = 0
        k = 0

        for cls in all_classes:
            axioms = chebi_ontology.get_axioms_for_iri(
                cls,
            )  # Get all axioms for the class
            cls_str = str(cls)
            is_deprecated = False
            has_smiles = False
            smiles_value = None

            for axiom in axioms:
                component = axiom.component  # The component of the axiom can e.g. be SubClassOf, AnnotationAssertion, etc.
                if (
                    type(component).__name__ == "AnnotationAssertion"
                ):  # Check if axiom is an annotation, e.g. a SMILES or deprecated tag
                    ann = component.ann  # Get the annotation
                    if hasattr(
                        ann,
                        "ap",
                    ):  # Check if annotation has an annotation property (ap). This can e.g. tell us if it is a SMILES string
                        prop_str = str(
                            ann.ap,
                        )  # Converts property IRI to string for easier comparison

                        if prop_str == deprecated_prop_str:  # Check if deprecated
                            is_deprecated = True
                            k += 1
                            if k % 1000 == 0:
                                print(f"Found {k} deprecated classes so far...")

                        if prop_str == smiles_prop_str:  # Check if SMILES property
                            has_smiles = True
                            smiles_value = str(ann.av)
                            i += 1
                            if i % 25000 == 0:
                                print(f"Found {i} classes with SMILES so far...")

            # Add to correct lists
            if is_deprecated:
                deprecated_classes.append(cls_str)
                continue  # Skip further checks for deprecated classes
            if has_smiles:
                classes_with_smiles.append(cls_str)
                # A class is a leaf if it has its own valid, non-wildcard SMILES, regardless of
                # whether it has subclasses. Classes with children that qualify are spliced out
                # of the hierarchy below so leaves stay terminal. Wildcard/R-group placeholder
                # SMILES (e.g. ChEBI's '*') are excluded.
                if not _has_wildcard(smiles_value):
                    leaf_classes_with_smiles.append(cls_str)

                    j += 1
                    if j % 10000 == 0:
                        print(f"Found {j} leaf classes with SMILES so far...")
            # Helper function to get all ancestors recursively

    # Sanity guard: a near-empty leaf set means leaf detection silently failed (e.g. SMILES
    # values were not read as expected). Catch it here, because the splice verification below
    # passes trivially when there are no leaves to misplace.
    if len(classes_with_smiles) > 0 and len(leaf_classes_with_smiles) < 0.5 * len(
        classes_with_smiles,
    ):
        raise RuntimeError(
            f"Only {len(leaf_classes_with_smiles)} leaf classes found among "
            f"{len(classes_with_smiles)} classes with SMILES — leaf detection likely failed "
            "(check how SMILES annotation values are being read).",
        )

    # Flatten the hierarchy so leaves are terminal: every class reconnects to its nearest
    # non-leaf ancestors, climbing past any leaf parents (handles chains of stacked leaves).
    # After this no leaf is anyone's parent, so leaves become siblings under the real categories.
    leaf_set = set(leaf_classes_with_smiles)
    print("Splicing leaves out of the hierarchy so they become terminal siblings...")
    flattened_subclass_map, flattened_parent_map = _splice_leaves_from_hierarchy(
        direct_parent_map,
        leaf_set,
    )
    _verify_splice(direct_parent_map, flattened_parent_map, leaf_set)

    # Overwrite the saved subclass map with the flattened version so every downstream artifact
    # (parent map, all-ancestors map, class->leaf map, graph) inherits the sibling structure.
    with open(subclass_map_file, "w") as f:
        json.dump(flattened_subclass_map, f)
    print(f"Saved flattened subclass map to {subclass_map_file}.")

    # Build leaf to all (non-leaf) ancestors map from the flattened hierarchy
    leaf_to_parents = {}
    i = 0
    ancestor_cache_from_parent_map = {}
    for i, leaf in enumerate(leaf_classes_with_smiles, start=1):
        all_ancestors = _get_all_ancestors_from_parent_map(
            leaf,
            flattened_parent_map,
            ancestor_cache_from_parent_map,
        )
        leaf_to_parents[leaf] = all_ancestors
        if i % 25000 == 0:
            print(f"Processed {i} leaf classes for ancestor mapping...")

    # Save leaf to ALL ancestors map
    with open(leaf_parents_map_file, "w") as f:
        json.dump(leaf_to_parents, f, indent=2)
    print(
        f"Saved leaf to all ancestors map to {leaf_parents_map_file} with {len(leaf_to_parents)} entries.",
    )

    # Summary
    print(f"\nDeprecated classes: {len(deprecated_classes)}")
    print(f"Classes WITH SMILES: {len(classes_with_smiles)}")
    print(f"Leaf classes with SMILES: {len(leaf_classes_with_smiles)}")

    return set(leaf_classes_with_smiles), set(deprecated_classes)


def save_filtered_owl(
    chebi_file,
    classes_with_smiles_to_remove,
    deprecated_classes_to_remove,
    output_file,
):

    # Check if output file already exists to avoid overwriting
    try:
        with open(output_file):
            print(
                f"Output file {output_file} already exists. Please remove it or change name before running this function.",
            )
            return
    except FileNotFoundError:
        pass  # File does not exist, proceed

    print(f"\nSaving filtered ontology to {output_file}...")

    tree = ET.parse(
        chebi_file,
    )  # Parsing the original OWL file to an ElementTree to be able to remove elements etc
    root = tree.getroot()

    classes_to_remove = classes_with_smiles_to_remove.union(
        deprecated_classes_to_remove,
    )

    # Define namespaces
    ns = {
        "owl": "http://www.w3.org/2002/07/owl#",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    }

    # Collect elements to remove
    elements_to_remove = []
    i = 0
    j = 0

    for elem in list(root):
        # Check if element is a Class
        if elem.tag == f"{{{ns['owl']}}}Class":
            class_iri = elem.attrib.get(f"{{{ns['rdf']}}}about")  # Get the class IRI

            # If class NOT in keep list, mark for removal
            if class_iri and class_iri in classes_to_remove:
                elements_to_remove.append(elem)
                i += 1
                if i % 1000 == 0:
                    print(f"Marked {i} classes for removal so far...")

    print(f"Total elements marked for removal: {len(elements_to_remove)}")

    # Remove all marked elements
    # Could combine marking elements and removing them -- check for later
    for j, elem in enumerate(elements_to_remove, start=1):
        root.remove(elem)
        if j % 5000 == 0:
            print(f"Removed {j} elements so far...")

    # Save filtered OWL
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"✓ Done! Saved {output_file}")


def _build_smiles_map_from_owl(owl_file, smiles_property):
    """OPTIMIZATION: Parse OWL XML directly to build SMILES map (much faster than ontology API)."""
    smiles_map = {}

    ns = {
        "owl": "http://www.w3.org/2002/07/owl#",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    }

    annotation_tag = f"{{{ns['owl']}}}AnnotationAssertion"
    annotation_prop_tag = f"{{{ns['owl']}}}annotationProperty"
    annotation_target_tag = f"{{{ns['owl']}}}annotationSubject"
    annotation_value_tag = f"{{{ns['owl']}}}annotationValue"
    rdf_resource = f"{{{ns['rdf']}}}resource"

    context = ET.iterparse(owl_file, events=("end",))

    for event, elem in context:
        if elem.tag == annotation_tag:
            prop_elem = elem.find(annotation_prop_tag)
            if prop_elem is not None:
                prop_resource = prop_elem.get(rdf_resource, "")
                if smiles_property in prop_resource:
                    target_elem = elem.find(annotation_target_tag)
                    value_elem = elem.find(annotation_value_tag)
                    if target_elem is not None and value_elem is not None:
                        class_iri = target_elem.get(rdf_resource, "")
                        smiles = value_elem.text
                        if class_iri and smiles:
                            smiles_map[class_iri] = smiles
        elem.clear()

    return smiles_map


def save_leaf_classes_with_smiles(
    leaf_classes,
    chebi_ontology,
    smiles_property,
    output_file,
    structural_classes,
    functional_classes,
    owl_file=None,
    parent_map_file=None,
):
    """Save leaf classes with SMILES to a CSV file.

    OPTIMIZATION: If owl_file and parent_map_file are provided, uses fast lookups instead of slow ontology API.
    """

    print(f"\nSaving {len(leaf_classes)} leaf classes with SMILES to {output_file}...")

    # OPTIMIZATION: Load precomputed maps instead of using slow ontology API
    smiles_map = {}
    parent_map = {}

    if owl_file and os.path.exists(owl_file):
        print("  Using fast OWL XML parsing for SMILES...")
        smiles_map = _build_smiles_map_from_owl(owl_file, smiles_property)
        print(f"  Loaded {len(smiles_map)} SMILES entries")

    if parent_map_file and os.path.exists(parent_map_file):
        print("  Loading parent map from JSON...")
        with open(parent_map_file) as f:
            parent_map = json.load(f)
        print(f"  Loaded {len(parent_map)} parent relationships")

    seen = set()
    rows = []

    for cls in leaf_classes:
        if cls in seen:
            print(f"Skipping duplicate class: {cls}")
            continue  # Skip duplicates
        seen.add(cls)

        # Try OWL map first, fall back to ontology API
        smiles = smiles_map.get(cls)
        if smiles is None and chebi_ontology:
            axioms = chebi_ontology.get_axioms_for_iri(cls)
            for axiom in axioms:
                component = axiom.component
                if type(component).__name__ == "AnnotationAssertion":
                    ann = component.ann
                    if hasattr(ann, "ap"):
                        prop_str = str(ann.ap)
                        if prop_str == f"<{smiles_property}>":
                            smiles = str(ann.av)
                            break

        # Try parent map first, fall back to ontology API
        classification = "neither"
        parents = set(parent_map.get(cls, []))

        if not parents and chebi_ontology:
            try:
                parents = {str(p) for p in chebi_ontology.get_superclasses(cls)}
            except Exception as e:  # noqa: BLE001
                print(f"⚠️ Could not get superclasses for {cls}: {e}")
                classification = "unknown"
                parents = set()

        if parents:
            if any(p in structural_classes for p in parents):
                classification = "structural"
            elif any(p in functional_classes for p in parents):
                classification = "functional"

        rows.append([cls, smiles, classification])

    # Write to CSV
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["IRI", "SMILES", "Classification"])
        writer.writerows(rows)

    length = len(rows)
    print(f"✓ Done! Saved {length} leaf classes with SMILES to {output_file}")


# Only used when deprecated classes were removed in a separate step. Should not be needed now. # Used in build_parent_map
def find_deprecated_classes(ontology, deprecated_property):
    """Find all deprecated classes in the ontology."""
    print("Finding deprecated classes...")
    all_classes = ontology.get_classes()
    deprecated_classes = set()
    i = 0

    for cls in all_classes:
        axioms = ontology.get_axioms_for_iri(cls)
        for axiom in axioms:
            component = axiom.component
            if type(component).__name__ == "AnnotationAssertion":
                ann = component.ann
                if hasattr(ann, "ap"):
                    prop_str = str(ann.ap)
                    if prop_str == f"<{deprecated_property}>":
                        deprecated_classes.add(cls)
                        i += 1
                        if i % 1000 == 0:
                            print(f"Found {i} deprecated classes so far...")
                        break  # No need to check more axioms for this class
    print(f"Total deprecated classes found: {len(deprecated_classes)}")
    return deprecated_classes


# Only used when deprecated classes were removed in a separate step. Should not be needed now.
def remove_classes_from_owl(input_file, classes_to_remove, output_file):
    """Remove given classes from OWL file and save a new version."""
    print(f"\nRemoving {len(classes_to_remove)} classes from {input_file}...")

    tree = ET.parse(input_file)
    root = tree.getroot()

    ns = {
        "owl": "http://www.w3.org/2002/07/owl#",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    }

    removed = 0
    for elem in list(root):
        if elem.tag == f"{{{ns['owl']}}}Class":
            iri = elem.attrib.get(f"{{{ns['rdf']}}}about")
            if iri and iri in classes_to_remove:
                root.remove(elem)
                removed += 1
                if removed % 5000 == 0:
                    print(f"Removed {removed} deprecated classes so far...")

    tree.write(output_file, encoding="utf-8", xml_declaration=True)
    print(f"✓ Done! Removed {removed} classes.")
    print(f"Saved to {output_file}")


# For building parent map.


def build_parent_map(
    ontology,
    output_json,
    deprecated_property,
    subclass_map_file=None,
    precomputed_deprecated_classes=None,
):
    """Build parent map excluding deprecated classes.

    Fast path: if subclass_map_file is available, derive child->parents from that JSON.
    """

    # Find deprecated classes first
    if precomputed_deprecated_classes is not None:
        deprecated_classes = {str(c) for c in precomputed_deprecated_classes}
        print(f"Using precomputed deprecated classes: {len(deprecated_classes)}")
    else:
        raw_deprecated = find_deprecated_classes(ontology, deprecated_property)
        deprecated_classes = {str(c) for c in raw_deprecated}

    print(f"Excluding {len(deprecated_classes)} deprecated classes")

    # Convert deprecated classes to strings for consistent comparison

    parent_map = {}

    if subclass_map_file and os.path.exists(subclass_map_file):
        print(
            f"Loading subclass map from {subclass_map_file} for fast parent-map build...",
        )
        with open(subclass_map_file) as f:
            subclass_map = json.load(f)

        parent_map_full = _invert_subclass_map(subclass_map)
        for cls_id, parents in parent_map_full.items():
            if cls_id in deprecated_classes:
                continue
            parent_map[cls_id] = [p for p in parents if p not in deprecated_classes]
        print(f"Built parent map from subclass map for {len(parent_map)} classes")
    else:
        all_classes = ontology.get_classes()
        string_ids = {cls: str(cls) for cls in all_classes}  # fast lookup table

        i = 0
        for cls in all_classes:
            if i % 1000 == 0:
                print(f"Processed {i} classes for parent map...")
            i += 1

            cls_id = string_ids[cls]

            # Skip if deprecated (now comparing strings)
            if cls_id in deprecated_classes:
                continue

            # Get parents, excluding deprecated ones
            direct_parents = ontology.get_superclasses(cls)
            parents = [
                string_ids[parent]
                for parent in direct_parents
                if string_ids[parent] not in deprecated_classes
            ]

            parent_map[cls_id] = parents

    with open(output_json, "w") as f:
        json.dump(parent_map, f, indent=2)

    # Print statistics
    root_classes = [cls for cls, parents in parent_map.items() if not parents]
    print(f"Saved {len(parent_map)} non-deprecated classes")
    print(f"  - {len(root_classes)} root classes (no parents)")
    print(f"  - {len(parent_map) - len(root_classes)} classes with parents")
    print(f"Saved parent map to {output_json}")


def get_name(root, iri, ns):
    """Return the rdfs:label for a class IRI."""
    for cls in root.findall("owl:Class", ns):
        about = cls.attrib.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about")
        if about == iri:
            label_elem = cls.find("rdfs:label", ns)
            if label_elem is not None:
                return label_elem.text.strip()
    return None


def shorten_parent_map(input_json, output_json):
    """Shorten the IRIs in the parent map by removing the common prefix."""
    prefix = "http://purl.obolibrary.org/obo/"

    with open(input_json) as f:
        parent_map = json.load(f)

    short_map = {}

    for iri, parents in parent_map.items():
        # Convert key
        short_key = iri.replace(prefix, "")

        # Convert parent IRIs
        short_parents = [p.replace(prefix, "") for p in parents]

        short_map[short_key] = short_parents

    with open(output_json, "w") as f:
        json.dump(short_map, f, indent=2)

    print(f"Saved shortened parent map → {output_json}")


def map_names_to_classes(chebi_ontology, output_json):
    ns = {
        "owl": "http://www.w3.org/2002/07/owl#",
        "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    }

    tree = ET.parse(chebi_ontology)
    root = tree.getroot()

    prefix = "http://purl.obolibrary.org/obo/"

    iri_to_name = {}
    i = 0

    for cls in root.findall("owl:Class", ns):
        iri = cls.attrib.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about")
        if iri is None:
            continue

        # keep only the last part, e.g. CHEBI_12345
        if iri.startswith(prefix):
            short_id = iri[len(prefix) :]
        else:
            short_id = iri  # fallback

        label_elem = cls.find("rdfs:label", ns)
        name = (
            label_elem.text.strip()
            if label_elem is not None and label_elem.text is not None
            else None
        )

        iri_to_name[short_id] = name

        i += 1
        if i % 25000 == 0:
            print(f"Processed {i} classes for name mapping...")

    with open(output_json, "w") as f:
        json.dump(iri_to_name, f, indent=2)

    print(f"Saved short-ID → name mapping to {output_json}")


if __name__ == "__main__":
    task = "build_parent_map"  # Options: "remove_leaves_with_smiles", "save_removed_leaf_classes", "build_parent_map", "map_names_to_classes"
    print(f"Selected task: {task}")
    # property IRI
    smiles_property = "https://w3id.org/chemrof/smiles_string"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"

    if (
        task == "remove_leaves_with_smiles"
    ):  # Remove leaf classes with SMILES and deprecated classes from OWL
        # (a new filtered OWL file will be saved)
        filtered_output_file = (
            "data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl"
        )

        chebi_file = "data/chebi.owl"
        subclass_map_file = "data/chebi_subclass_map.json"
        leaf_parents_map_file = "data/removed_leaf_classes_to_ALL_parents_map.json"

        use_found_leaf_classes = True  # Set to True to use previously found leaf classes with SMILES from CSV file
        removed_leaf_classes_file = "data/removed_leaf_classes_with_smiles.csv"  # Only needed if use_found_leaf_classes is True

        if os.path.exists(filtered_output_file):
            print(
                f"Output file {filtered_output_file} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.",
            )
        else:
            print(
                "Running code to remove leaf classes with SMILES and deprecated classes...",
            )
            chebi_ontology = load_chebi()
            classes_with_smiles, deprecated_classes = (
                find_leaf_classes_with_smiles_and_deprecated(
                    chebi_ontology,
                    smiles_property,
                    deprecated_property,
                    subclass_map_file,
                    leaf_parents_map_file,
                    use_found_leaf_classes,
                    removed_leaf_classes_file,
                )
            )

            # Comment out the next two lines if you do not want to save the filtered OWL in a new file
            save_filtered_owl(
                chebi_file,
                classes_with_smiles,
                deprecated_classes,
                filtered_output_file,
            )
            filtered_ontology = load_ontology(
                filtered_output_file,
            )  # just to confirm it loads

    elif task == "save_removed_leaf_classes":
        print("Running code to save removed leaf classes with SMILES...")

        output_file = "data/removed_leaf_classes_with_smiles_new.csv"
        subclass_map_file = "data/chebi_subclass_map.json"

        structural_ontology = load_ontology(
            "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl",
        )
        functional_ontology = load_ontology(
            "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl",
        )

        # Extract sets of class IRIs for quick membership checking
        structural_classes = set(structural_ontology.get_classes())
        functional_classes = set(functional_ontology.get_classes())

        print(f"Structural classes loaded: {len(structural_classes)}")
        print(f"Functional classes loaded: {len(functional_classes)}")

        chebi_ontology = load_chebi()
        classes_with_smiles, _, _ = find_leaf_classes_with_smiles_and_deprecated(
            chebi_ontology,
            smiles_property,
            deprecated_property,
            subclass_map_file,
        )

        # Save CSV with classification
        save_leaf_classes_with_smiles(
            classes_with_smiles,
            chebi_ontology,
            smiles_property,
            output_file,
            structural_classes,
            functional_classes,
        )

        df = pd.read_csv(output_file)
        print(len(df), "rows total")
        print(df["IRI"].nunique(), "unique IRIs")
        print(df[df["SMILES"].isna()])

    elif task == "build_parent_map":
        print("Running code to build parent map excluding deprecated classes...")
        ontology = load_chebi()
        output_json = "data/chebi_parent_map.json"
        # shortened_output_json = "data/chebi_parent_map_shortened_id.json"
        if os.path.exists(output_json):
            print(
                f"Output file {output_json} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.",
            )
        else:
            build_parent_map(
                ontology,
                output_json,
                deprecated_property,
            )  # Deprecated classes are removed inside the function

    elif task == "map_names_to_classes":
        print("Running code to map names to classes...")
        chebi_ontology = "data/chebi.owl"
        output_json = "data/chebi_id_to_name_map.json"
        if os.path.exists(output_json):
            print(
                f"Output file {output_json} already exists. Are you sure you want to overwrite it? If so, please remove it before running this script.",
            )
        else:
            map_names_to_classes(chebi_ontology, output_json)

    else:
        print("No valid task selected. Please choose a valid task.")
