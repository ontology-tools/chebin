"""Count classes removed from a given class in the structural split.

The total number of classes in the structural split is also counted.
OBS: make sure to only count a removed class once.
"""

import json
import logging
from collections import defaultdict

import pandas as pd

logger = logging.getLogger(__name__)
_MISSING_ROLE_WARNING_LIMIT = 5
_missing_role_warning_count = 0
_missing_role_warning_suppressed = 0

_structural_leaf_ids_cache = {}


def get_structural_leaf_ids(removed_classes_csv: str) -> set[str]:
    """
    Return the set of leaf class IRIs correctly classified as 'structural'.
    Leaves classified as 'functional' or 'neither' are ChEBI ontology labeling
    mistakes (a genuine leaf class should only exist under structural), so they
    are excluded everywhere a background/study leaf set is built.
    """
    if removed_classes_csv not in _structural_leaf_ids_cache:
        removed_classes = pd.read_csv(removed_classes_csv)
        structural = removed_classes[removed_classes["Classification"] == "structural"]
        _structural_leaf_ids_cache[removed_classes_csv] = set(structural["IRI"])
    return _structural_leaf_ids_cache[removed_classes_csv]


def count_removed_leaves(removed_classes_csv: str) -> int:
    return len(get_structural_leaf_ids(removed_classes_csv))


def build_class_to_leaf_map(
    leaf_to_ancestors_file: str,
    class_to_leaf_output_file: str,
) -> None:
    """Build a JSON map from each class IRI to ALL its leaf descendants using an existing leaf-to-ancestors map."""

    print(f"Loading leaf to ancestors map from {leaf_to_ancestors_file}...")
    with open(leaf_to_ancestors_file) as f:
        leaf_to_ancestors = json.load(f)
    print(f"Loaded leaf to ancestors map with {len(leaf_to_ancestors)} leaf classes.")

    class_to_leaves: dict[str, set[str]] = defaultdict(set)
    print("Building class to leaf descendants map...")

    # For each leaf, add it to all its ancestors
    for leaf, ancestors in leaf_to_ancestors.items():
        for ancestor in ancestors:
            class_to_leaves[ancestor].add(leaf)

    # Convert sets to lists for JSON serialization
    class_to_leaf_json = {cls: list(leaves) for cls, leaves in class_to_leaves.items()}

    # Save to JSON
    with open(class_to_leaf_output_file, "w") as f:
        json.dump(class_to_leaf_json, f, indent=2)

    print(
        f"Saved class to leaf descendants map with {len(class_to_leaf_json)} classes to {class_to_leaf_output_file}.",
    )


def count_removed_classes_for_class(
    class_iri: str,
    class_to_leaf_map: dict[str, list[str]],
    classification: str,
    class_to_all_roles_map: dict[str, list[str]],
    roles_to_leaves_map: dict[str, list[str]],
    structural_leaf_ids: set[str] | None = None,
) -> tuple[set[str], int]:
    """
    Count how many leaf classes are associated with the given class_iri.

    Parameters:
        class_iri: The class to check
        class_to_leaf_map: Maps structural classes to their leaf descendants
        classification: "structural", "functional", or "full"
        class_to_all_roles_map: Maps classes to all roles (used for functional classification)
        roles_to_leaves_map: Maps role classes to their associated leaves
        structural_leaf_ids: If given, leaf descendants are restricted to this
            set (the genuine 'structural'-classified leaves), excluding any
            leaf classes mislabeled with another Classification in ChEBI

    Returns:
        Tuple of (leaves set, leaf count).

    Raises:
        ValueError: If class_iri not in class_to_leaf_map for structural/full classification
    """

    if classification not in ["structural", "functional", "full"]:
        raise ValueError(
            f"Classification must be 'structural', 'functional', or 'full', got: {classification}",
        )

    leaves: set[str] = set()

    if classification in ["structural", "full"]:
        if str(class_iri) not in class_to_leaf_map:
            raise ValueError(
                f"Class {class_iri} not found in class_to_leaf_map. "
                "Check that the map file was loaded correctly and the class IRI is valid.",
            )

        leaves_structure = class_to_leaf_map[str(class_iri)]
        if len(leaves_structure) == 0:
            raise ValueError(
                f"Class {class_iri} has no leaf descendants in ontology. "
                "This should not happen if the map was built correctly.",
            )

        leaves.update(leaves_structure)

        if structural_leaf_ids is not None:
            leaves &= structural_leaf_ids

    else:
        raise ValueError("Functional classification is not yet implemented.")

    n_leaves = len(leaves)

    return leaves, n_leaves


def count_removed_classes_for_roles(
    class_iri: str,
    leaf_descendants_map: dict[str, list[str]],
    classification: str,
    roles_to_leaves_map: dict[str, list[str]],
) -> tuple[set[str], int]:
    """
    Count leaf classes associated with the given role IRI.

    Parameters:
        class_iri: The role IRI to check
        leaf_descendants_map: Unused (kept for API compatibility)
        classification: "functional" or "full"
        roles_to_leaves_map: Maps role classes to their associated leaves

    Returns:
        Tuple of (leaves set, leaf count).
    """
    if classification not in ["functional", "full"]:
        print(f"Classification {classification} is not supported for counting roles.")
        return set(), 0

    leaves: set[str] = set()
    leaves.update(roles_to_leaves_map.get(class_iri, []))

    if not leaves:
        print(f"⚠️ Role {class_iri} has no associated leaves in roles_to_leaves_map.")

    n_leaves = len(leaves)

    return leaves, n_leaves


if __name__ == "__main__":
    """ Select task to perform. To calculate the removed leaf classes for a given class,
    the file with the class to leaf descendants map must be created first.
    This is done in "build_class_to_leaf_map" task. """

    task = "count_removed_classes_for_class"
    # Options: "count_total_removed_leaves" "count_removed_classes_for_class" "build_class_to_leaf_map" "enrichment_analysis_plain"

    # Variables used in "count_removed_classes_for_class":
    # - classification (only if check_leaf_classes is True), check_leaf_classes, class_iri

    # Variables used in "count_total_removed_leaves":
    # - classification

    # Variables used in "enrichment_analysis_plain":
    # - classification, class_iri, (check_leaf_classes)

    # Relevant for task "count_removed_classes_for_class"
    classification = "structural"  # "functional" or "structural" or "full"
    check_leaf_classes = True  # Checks that the found the leaf classes are of the exppected type (Functional or Structural)
    class_iri = "http://purl.obolibrary.org/obo/CHEBI_83822"

    # Files
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
    map_file = "data/class_to_leaf_descendants_map.json"
    leaf_to_ancestors_file = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_output_file = "data/class_to_leaf_descendants_map.json"
    class_to_all_roles_json = "data/class_to_all_roles_map.json"
    roles_to_leaves_map_json = "data/roles_to_leaves_map.json"

    with open(class_to_all_roles_json) as f:
        class_to_all_roles_map = json.load(f)
    with open(roles_to_leaves_map_json) as f:
        roles_to_leaves_map = json.load(f)

    if task == "count_total_removed_leaves":
        n_removed_leaves = count_removed_leaves(removed_leaves_csv)
        print(
            f"Total number of removed {classification} leaf classes: {n_removed_leaves}",
        )

        # Output
        # - Total number of removed structural leaf classes: 184 436
        # - Total number of removed functional leaf classes: 41

    elif task == "count_removed_classes_for_class":
        with open(map_file) as f:
            class_to_leaf_map = json.load(f)
        subclasses, n_subclasses = count_removed_classes_for_class(
            class_iri,
            class_to_leaf_map,
            classification,
            class_to_all_roles_map,
            roles_to_leaves_map,
        )
        print(
            f"Class {class_iri} has {n_subclasses} removed subclasses in the ontology.",
        )

        if subclasses and n_subclasses < 50:
            for sub in subclasses:
                print(f" - {sub}")

    elif task == "build_class_to_leaf_map":
        build_class_to_leaf_map(leaf_to_ancestors_file, class_to_leaf_output_file)

    elif task == "enrichment_analysis_plain":  # to be removed
        _, n_subclasses = count_removed_classes_for_class(
            class_iri,
            map_file,
            classification,
            class_to_all_roles_map,
            roles_to_leaves_map,
        )
        n_tot_removed_leaves = count_removed_leaves(removed_leaves_csv)

        print(
            f"Calculated for {classification} ontology for {class_iri} (make sure classification is correct):",
        )
        prob_of_success = n_subclasses / n_tot_removed_leaves

        print(
            f" The probability of success value for class {class_iri} is: {prob_of_success:.4g} ({n_subclasses} / {n_tot_removed_leaves})",
        )

    else:
        print("No valid task selected.")
