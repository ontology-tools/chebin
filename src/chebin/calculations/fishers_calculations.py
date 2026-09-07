"""Fisher-based enrichment calculations and helpers."""

import json
import time

import pandas as pd
from scipy.stats import fisher_exact

from chebin.calculations.chebi_ids import CHEBI_IRI_PREFIX, to_chebi_iri
from chebin.calculations.data_files import load_data_files
from chebin.calculations.log_utils import describe
from chebin.calculations.multiple_test_corrections import (
    benjamini_hochberg_fdr_correction,
    bonferroni_correction,
)
from chebin.calculations.pre_fishers_calculations import (
    count_removed_classes_for_class,
    count_removed_classes_for_roles,
    count_removed_leaves,
    get_structural_leaf_ids,
)
from chebin.calculations.smiles_lookup import smiles_list_to_studyset
from chebin.calculations.visualitations_and_pruning import (
    create_graph_with_roles_and_structures,
    high_p_value_branch_pruner,
    id_to_name,
    linear_branch_collapser_pruner_remove_less,
    root_children_pruner,
    zero_degree_pruner,
)

# Contingency table:
#   a = study compounds annotated to this class
#   b = study compounds not annotated to this class
#   c = background compounds in this class (excluding the study set)
#   d = background compounds not in this class (excluding the study set)


def calculate_p_value(
    n_ss_annotated: int,
    n_ss_leaves: int,
    n_bg_annotated: int,
    n_bg_leaves: int,
) -> tuple[float | None, float | None]:
    """
    Calculate Fisher's exact test p-value for enrichment.

    Constructs a 2x2 contingency table from enrichment counts and returns
    the odds ratio and p-value using Fisher's exact test (greater-tail).

    Contingency table:
        - a: items in study set annotated to class
        - b: items in study set NOT annotated to class
        - c: items in background NOT in study set, annotated to class
        - d: items in background NOT in study set, NOT annotated to class

    Args:
        n_ss_annotated: Count of study set items annotated to this class.
        n_ss_leaves: Total count of items in study set.
        n_bg_annotated: Count of background items annotated to this class.
        n_bg_leaves: Total count of items in background.

    Returns:
        tuple: (odds_ratio, p_value) if valid contingency table, else (None, None).
            Odds ratio and p-value for greater-tail Fisher's exact test.
    """
    a = n_ss_annotated
    b = n_ss_leaves - n_ss_annotated
    c = n_bg_annotated - n_ss_annotated
    d = n_bg_leaves - n_bg_annotated - b

    # Validate contingency table values are non-negative
    if a < 0 or b < 0 or c < 0 or d < 0:
        print(
            f"Warning: Invalid contingency table with negative values: a={a}, b={b}, c={c}, d={d}",
        )
        print(
            f"  (n_ss_annotated={n_ss_annotated}, n_ss_leaves={n_ss_leaves}, n_bg_annotated={n_bg_annotated}, n_bg_leaves={n_bg_leaves})",
        )
        return None, None

    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    # odds, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
    # ‘greater’: odds ratio of the underlying population is greater than one
    return odds, p


def normalize_id(raw_id: str) -> str:
    """
    Normalize ChEBI identifiers to full IRIs.

    Converts every ChEBI ID format users write (CHEBI:12345, chebi:12345,
    CHEBI ID: 12345, CHEBI_12345, the bare number, or a full IRI) to the standard
    http://purl.obolibrary.org/obo/CHEBI_XXXXX IRI format. See
    :mod:`chebin.calculations.chebi_ids` for the forms recognised.

    Args:
        raw_id: ChEBI identifier in any format (CHEBI:12345, numeric ID, or IRI).

    Returns:
        str: Normalized IRI in http://purl.obolibrary.org/obo/CHEBI_XXXXX format.
    """
    iri = to_chebi_iri(raw_id)
    if iri is not None:
        return iri
    # Not a ChEBI ID in any recognised form: keep the historical behaviour of
    # passing IRIs through untouched and putting anything else in the OBO
    # namespace, so unrecognised input still fails as a non-matching node rather
    # than as an exception here.
    value = raw_id.strip().replace('"', "")
    if value.startswith(("http://", "https://")):
        return value
    return f"{CHEBI_IRI_PREFIX}{value.replace(':', '_')}"


def get_leaves(studyset_list, leaves_csv, class_to_leaf_map, structural_leaf_ids=None):
    """Get leaf descendants for the input classes, or the class itself if it's already a leaf.

    structural_leaf_ids: if given, leaves mislabeled with a Classification other
    than 'structural' (a ChEBI ontology error) are excluded.
    """
    studyset_leaves = set()

    leaves_df = pd.read_csv(leaves_csv)
    all_leaf_ids = set(leaves_df["IRI"].values)

    for cls in studyset_list:
        print(f"Processing class {cls}...")
        if cls in all_leaf_ids:
            if structural_leaf_ids is not None and cls not in structural_leaf_ids:
                print(
                    f"Excluding class {cls}: not classified as 'structural' in ChEBI (likely a mislabeled leaf).",
                )
                continue
            # add to list of studyset leaves
            print(f"Class {cls} is already leaf.")
            studyset_leaves.add(cls)
        else:
            # get leaf descendants from map and add them to studyset leaves
            leaf_descendants = set(class_to_leaf_map.get(cls, []))
            if structural_leaf_ids is not None:
                leaf_descendants &= structural_leaf_ids
            print(
                f"Class {cls} is not a leaf, adding its {len(leaf_descendants)} leaf descendants.",
            )
            studyset_leaves.update(leaf_descendants)

    return list(studyset_leaves)


def get_ancestors_for_inputs(
    studyset_leaves: list[str],
    leaf_to_all_parents_map_json: str,
) -> list[str]:
    """
    Extract all ancestors (parents at all levels) of the given leaf classes.

    Args:
        studyset_leaves: Iterable of leaf class IRIs.
        leaf_to_all_parents_map_json: Path to JSON file mapping leaf IRIs to
            lists of all ancestor (parent) IRIs.

    Returns:
        list: All unique ancestor class IRIs reachable from studyset_leaves.
    """
    with open(leaf_to_all_parents_map_json) as f:
        leaf_to_all_parents_map = json.load(f)

    studyset_ancestors: set[str] = set()
    for leaf in studyset_leaves:
        parents = leaf_to_all_parents_map.get(leaf, [])
        studyset_ancestors.update(parents)

    return list(studyset_ancestors)


def get_n_ss_annotated(
    studyset_leaves,
    class_to_check,
    class_to_leaf_map,
    classification,
    class_to_all_roles_map,
    roles_to_leaves_map,
):
    """
    n_ss_annotated = number of input classes that are leaf descendants of the given class.

    Parameters:
        studyset_leaves: list of class IDs that were provided by the user (the study set)
        class_to_check: the ontology class for which we want n_ss_annotated
        class_to_leaf_map: JSON mapping each class to all its leaf descendants
        classification: "structural", "functional", or "full"
        class_to_all_roles_map: Maps classes to all roles (used for functional classification)
        roles_to_leaves_map: Maps role classes to their associated leaves

    Returns:
        int: Count of study set leaves that are descendants of class_to_check

    Raises:
        ValueError: If classification is not supported or class not found
    """

    leaves = set()

    if classification in ["structural", "full"]:
        # descendants of the class we are calculating enrichment for
        if class_to_check not in class_to_leaf_map:
            raise ValueError(
                f"Class {class_to_check} not found in class_to_leaf_map. "
                "Ensure map file is loaded and class IRI is valid.",
            )

        leaf_descendants = set(class_to_leaf_map.get(class_to_check, []))
        leaves.update(leaf_descendants)

        # count how many study classes appear in those leaf descendants
        n_ss_annotated = len(leaves.intersection(set(studyset_leaves)))
        return n_ss_annotated

    # if classification in ["functional", "full"]:
    #     # Get all roles (direct + inherited from ancestors + role ancestors)
    #     all_roles = class_to_all_roles_map.get(class_to_check, [])
    #     for role in all_roles:
    #         leaves.update(roles_to_leaves_map.get(role, []))

    else:
        raise ValueError(
            f"Classification '{classification}' is not supported. "
            "Use 'structural', 'functional', or 'full'.",
        )


def get_n_ss_annotated_for_roles(
    studyset_leaves,
    class_to_check,
    class_to_all_roles_map,
    roles_to_leaves_map,
):
    """
    Count study set items annotated to a role/functional class.

    Args:
        studyset_leaves: Set/list of study set leaf class IRIs.
        class_to_check: Role/functional class IRI to check annotation for.
        class_to_all_roles_map: Map (unused in this function, kept for API consistency).
        roles_to_leaves_map: Mapping from role IRIs to lists of leaf IRIs.

    Returns:
        int: Number of study set items with this role annotation.
    """
    leaves = set()
    leaves.update(roles_to_leaves_map.get(class_to_check, []))
    n_ss_annotated = len(leaves.intersection(set(studyset_leaves)))
    return n_ss_annotated


def get_enrichment_values(
    removed_leaves_csv,
    classification,
    studyset_leaves,
    studyset_ancestors,
    class_to_leaf_map,
    class_to_all_roles_map,
    roles_to_leaves_map,
    studyset_ancestors_roles,
    structural_leaf_ids=None,
):
    """
    Calculate Fisher's exact p-values for enrichment of classes in study set.

    Iterates over all candidate classes (structural ancestors or functional roles)
    and computes the contingency table for each, returning p-values and odds ratios.

    Args:
        removed_leaves_csv: Path to CSV with all leaf classes and their counts.
        classification: Classification type: "structural" or "functional".
        studyset_leaves: List of leaf class IRIs in study set.
        studyset_ancestors: List of all ancestor class IRIs reachable from study set.
        class_to_leaf_map: Mapping from class IRIs to their leaf descendants.
        class_to_all_roles_map: Mapping from structural classes to functional roles.
        roles_to_leaves_map: Mapping from role IRIs to leaf descendants.
        studyset_ancestors_roles: List of role IRIs reachable from study set
            (only used if classification="functional").
        structural_leaf_ids: Optional set of valid structural leaf IRIs
            (filtering parameter).

    Returns:
        dict: Mapping from class/role IRIs to (odds_ratio, p_value) tuples.
    """

    # n_bg_leaves and n_ss_leaves will be the same for all classes
    n_bg_leaves = count_removed_leaves(removed_leaves_csv)
    n_ss_leaves = len(studyset_leaves)

    results = {}  # dictionary to hold results

    if classification in ["structural", "full"]:
        # Calculate enrichment for structural ancestors
        for class_to_check in studyset_ancestors:
            # print(f"Calculating enrichment for class {class_to_check}...")

            _, n_bg_annotated = count_removed_classes_for_class(
                class_to_check,
                class_to_leaf_map,
                classification,
                class_to_all_roles_map,
                roles_to_leaves_map,
                structural_leaf_ids,
            )
            n_ss_annotated = get_n_ss_annotated(
                studyset_leaves,
                class_to_check,
                class_to_leaf_map,
                classification,
                class_to_all_roles_map,
                roles_to_leaves_map,
            )

            odds, p_value = calculate_p_value(
                n_ss_annotated,
                n_ss_leaves,
                n_bg_annotated,
                n_bg_leaves,
            )

            # Skip if calculation failed due to invalid counts
            if p_value is None:
                print(
                    f"Skipping class {id_to_name(class_to_check)} due to invalid contingency table",
                )
                continue

            results[class_to_check] = {
                "class": id_to_name(class_to_check),
                # "class_id": strip_prefix(class_to_check),
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": odds,
                "p_value": p_value,
            }

        # if classification == "functional" or classification == "full":
        # # Update studyset_ancestors_roles with the roles associated (direct + inherited from ancestors) with the current class being checked
        # # These roles will be added to the graph
        #     studyset_ancestors_roles.update(class_to_all_roles_map.get(class_to_check, []))

    # Calculate enrichment for role classes
    if classification in ["functional", "full"] and studyset_ancestors_roles:
        print(f"Calculating enrichment for {len(studyset_ancestors_roles)} roles...")

        for role_to_check in studyset_ancestors_roles:
            _, n_bg_annotated = count_removed_classes_for_roles(
                role_to_check,
                class_to_leaf_map,
                classification,
                roles_to_leaves_map,
            )
            n_ss_annotated = get_n_ss_annotated_for_roles(
                studyset_leaves,
                role_to_check,
                class_to_all_roles_map,
                roles_to_leaves_map,
            )

            odds, p_value = calculate_p_value(
                n_ss_annotated,
                n_ss_leaves,
                n_bg_annotated,
                n_bg_leaves,
            )

            # Skip if calculation failed due to invalid counts
            if p_value is None:
                print(
                    f"Skipping role {id_to_name(role_to_check)} due to invalid contingency table",
                )
                continue

            results[role_to_check] = {
                "class": id_to_name(role_to_check),
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": odds,
                "p_value": p_value,
            }

    return results


def print_enrichment_results(enrichment_results):
    # Include corrected p-values if they exist
    # print(f"{'Class:':45} {'p-value:':15} {'p-value (corrected):':20} {'n_ss_annotated':20} {'n_bg_annotated':20}" )
    print(f"{'Class:':45} {'raw p-value:':15} {'p-value (corrected):':20}")
    print("-" * 200)

    for r in enrichment_results.values():
        corrected_p = r.get("p_value_corrected", None)
        corrected_p_str = f"{corrected_p:.4e}" if corrected_p is not None else "N/A"
        # print(f"{r['class']:45} {r['p_value']:.4e}      {corrected_p_str:20} {r['n_ss_annotated']:20} {r['n_bg_annotated']:20}")
        print(f"{r['class']:45} {r['p_value']:.4e}      {corrected_p_str:20}")


def write_enrichment_csv(enrichment_results: dict, csv_path: str) -> None:
    df = pd.DataFrame(list(enrichment_results.values()))
    columns = [
        col for col in ["class", "p_value", "p_value_corrected"] if col in df.columns
    ]
    df[columns].to_csv(csv_path, index=False)


# Graphing and pruning strategies.


def _load_data_files():
    """
    Load common data files needed for enrichment analysis.

    Returns:
        tuple: (class_to_leaf_map, class_to_all_roles_map, roles_to_leaves_map,
                removed_leaves_csv, leaf_to_ancestors_map_file, parent_map_file)
    """
    return load_data_files()


def run_enrichment_analysis(
    studyset_list: list[str],
    bonferroni_correct: bool = False,
    benjamini_hochberg_correct: bool = True,
    root_children_prune: bool = False,
    levels: int = 2,
    linear_branch_prune: bool = False,
    n: int = 2,
    high_p_value_prune: bool = False,
    p_value_threshold: float = 0.05,
    zero_degree_prune: bool = False,
    classification: str = "structural",
    print_results: bool = False,
    csv_output_path: str | None = None,
) -> tuple[dict, object]:
    """
    Run enrichment analysis with optional pruning and multiple test correction.

    Args:
        studyset_list: List of ChEBI class IDs (IRIs) to analyse.
        bonferroni_correct: Apply Bonferroni correction to p-values (bool, the standard is False).
        benjamini_hochberg_correct: Apply Benjamini-Hochberg FDR correction to p-values (bool, the standard is True).
        root_children_prune: Apply root children pruning before enrichment (bool).
        levels: Number of levels to prune from root (int, default 2). Only used if root_children_prune is True.
        linear_branch_prune: Apply linear branch pruning before enrichment (bool).
        n: Keep every n-th node in linear branches (int, default 2). Only used if linear_branch_prune is True.
        high_p_value_prune: Apply high p-value pruning after enrichment (bool).
        p_value_threshold: Threshold for high p-value pruning (float, default 0.05). Only used if high_p_value_prune is True.
        zero_degree_prune: Apply zero-degree pruning after enrichment, removing nodes with no connections (bool).
        classification: Classification type for enrichment analysis ("structural", "functional", or "full").
        print_results: Print a results table to stdout before returning (bool, default False).
        csv_output_path: If given, write the enrichment results to this CSV path (default None).

    Returns:
        tuple: (results_dict, pruned_graph) where:
            - results_dict contains keys: "study_set" (item names), "removed_nodes"
              (pruned node names), "enrichment_results" (node -> p-value mapping)
            - pruned_graph is a networkx graph after all pruning operations
    """

    # Provide a warning if both Bonferroni and Benjamini-Hochberg corrections are requested
    if bonferroni_correct and benjamini_hochberg_correct:
        print(
            "Warning: Both Bonferroni and Benjamini-Hochberg corrections requested. "
            "Only one correction method should be applied at a time. "
            "Proceeding with Benjamini-Hochberg correction.",
        )
        bonferroni_correct = False

    pruning_before_enrichment = root_children_prune or linear_branch_prune

    # Load common data files
    (
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        removed_leaves_csv,
        leaf_to_ancestors_map_file,
        parent_map_file,
    ) = _load_data_files()

    structural_leaf_ids = get_structural_leaf_ids(removed_leaves_csv)

    # Normalize study set IDs to full IRIs
    studyset_list = [normalize_id(cls) for cls in studyset_list]

    studyset_leaves = get_leaves(
        studyset_list,
        removed_leaves_csv,
        class_to_leaf_map,
        structural_leaf_ids,
    )
    # print(f"Study set leaves: {studyset_leaves}")

    studyset_ancestors_all = get_ancestors_for_inputs(
        studyset_leaves,
        leaf_to_ancestors_map_file,
    )

    # print(f"Study set ancestors: {studyset_ancestors_all}")
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

    if classification in ["functional", "full"]:
        # Collect roles associated with structural ancestors AND leaves
        studyset_ancestors_roles = set()

        # Add roles from leaves
        for leaf in studyset_leaves:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(leaf, []))

        # Add roles from structural ancestors
        for class_to_check in studyset_ancestors_all:
            studyset_ancestors_roles.update(
                class_to_all_roles_map.get(class_to_check, []),
            )

        # Also include all ancestors of these roles (they will appear in the graph)
        with open(parent_map_file) as f:
            parent_map = json.load(f)

        roles_with_ancestors = set(studyset_ancestors_roles)
        to_process = list(studyset_ancestors_roles)

        while to_process:
            role = to_process.pop(0)
            for parent in parent_map.get(role, []):
                if parent not in roles_with_ancestors:
                    roles_with_ancestors.add(parent)
                    to_process.append(parent)

        studyset_ancestors_roles = roles_with_ancestors
        print(f"Number of roles (including ancestors): {len(studyset_ancestors_roles)}")
    else:
        studyset_ancestors_roles = set()

    G = create_graph_with_roles_and_structures(
        studyset_leaves,
        studyset_ancestors_all,
        studyset_ancestors_roles,
        parent_map_file,
        class_to_all_roles_map,
        classification,
    )

    pruned_G = G.copy()

    all_removed_nodes = set()

    if pruning_before_enrichment:
        if root_children_prune:
            print(f"studyset_leaves: {describe(studyset_leaves)}")
            print(f"Root children pruner activated, pruning {levels} levels from root")
            time_start_total = time.time()
            pruned_G, removed_nodes, _execution_count = root_children_pruner(
                pruned_G,
                levels,
                allow_re_execution=False,
                execution_count=0,
            )
            time_end_total = time.time()
            print(
                f"Total time for root children pruning: {time_end_total - time_start_total} seconds",
            )
            # print(f"Removed nodes by root children pruner: {describe(removed_nodes)}")
            all_removed_nodes.update(removed_nodes)

        if linear_branch_prune:
            print(
                f"Linear branch pruner activated, keeping only every {n}-th node in linear branches",
            )

            pruned_G, removed_nodes = linear_branch_collapser_pruner_remove_less(
                pruned_G,
                n,
            )
            # print(f"Removed nodes by linear branch pruner: {describe(removed_nodes)}")
            all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from studyset_ancestors_all
        studyset_ancestors = [
            cls for cls in studyset_ancestors_all if cls not in all_removed_nodes
        ]
        print(
            f"Number of study set ancestors after before-enrichment pruning: {len(studyset_ancestors)}",
        )

    else:
        studyset_ancestors = studyset_ancestors_all

    enrichment_results = get_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles,
        structural_leaf_ids,
    )

    # print("Enrichment results:")
    # print_enrichment_results(enrichment_results)

    if bonferroni_correct:
        print("Applying Bonferroni correction to p-values...")
        enrichment_results, _correction_map = bonferroni_correction(enrichment_results)
        # print("Enrichment results after Bonferroni correction:")
        # print_enrichment_results(enrichment_results)
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction to p-values...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
        # print("Enrichment results after Benjamini-Hochberg correction:")
        # print_enrichment_results(enrichment_results)

    if high_p_value_prune:  # Uses corrected p-values if correction was applied
        print(
            f"High p-value pruner activated, pruning nodes with p-value above {p_value_threshold}",
        )

        time_start_total = time.time()
        pruned_G, removed_nodes = high_p_value_branch_pruner(
            pruned_G,
            enrichment_results,
            p_value_threshold,
        )
        time_end_total = time.time()
        print(
            f"Total time for high p-value pruning: {time_end_total - time_start_total} seconds",
        )
        # print(f"Removed nodes by high p-value pruner: {describe(removed_nodes)}")
        print(f"Number of pruned nodes: {len(removed_nodes)}")

        # update all_removed_nodes
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

        # print("Final enrichment results after high p-value pruning:")
        # print_enrichment_results(enrichment_results)

    if zero_degree_prune:
        print("Applying zero-degree pruner to remove nodes with zero degree...")

        time_start_total = time.time()
        pruned_G, removed_nodes = zero_degree_pruner(pruned_G)
        time_end_total = time.time()
        print(
            f"Total time for zero-degree pruning: {time_end_total - time_start_total} seconds",
        )
        # print(f"Removed nodes by zero-degree pruner: {describe(removed_nodes)}")
        print(f"Number of pruned nodes: {len(removed_nodes)}")

        # update all_removed_nodes
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

        # print("Final enrichment results after zero-degree pruning:")
        # print_enrichment_results(enrichment_results)

    # print("Final enrichment results:")
    # print_enrichment_results(enrichment_results)

    print(f"Number of removed nodes in total: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {
            id_to_name(cls): vals for cls, vals in enrichment_results.items()
        },
    }
    if print_results:
        print_enrichment_results(results["enrichment_results"])
    if csv_output_path:
        write_enrichment_csv(results["enrichment_results"], csv_output_path)
    return results, pruned_G


def run_enrichment_analysis_from_smiles(
    smiles_list: list[str],
    use_parents: bool = False,
    **kwargs,
) -> tuple[dict, object, dict]:
    """Run run_enrichment_analysis on a list of SMILES strings instead of ChEBI IDs.

    Each SMILES is resolved to ChEBI ID(s) via chebin.calculations.smiles_lookup
    before delegating to run_enrichment_analysis; see that function for the
    remaining arguments (passed through as **kwargs) and return value.

    Args:
        smiles_list: List of SMILES strings to analyse.
        use_parents: If a SMILES can't be resolved to a direct ChEBI ID, fall back
            to a remote classification call and use its direct parent ChEBI IDs.

    Returns:
        tuple: (results_dict, pruned_graph, smiles_diagnostics) where
            smiles_diagnostics is {"unresolved_smiles": [...], "ambiguous_matches": [...]}.
    """
    studyset_list, unresolved_smiles, ambiguous_matches = smiles_list_to_studyset(
        smiles_list,
        use_parents=use_parents,
    )
    results, pruned_G = run_enrichment_analysis(studyset_list, **kwargs)
    return (
        results,
        pruned_G,
        {
            "unresolved_smiles": unresolved_smiles,
            "ambiguous_matches": ambiguous_matches,
        },
    )


####################################
# Combine pruning strategies
####################################

# Plain Enrichment Pruning Strategy: For the pre-loop phase this strategy applies the High Value Branch Pruner (0.05),
# the Linear Branch Collapser Pruner, and the Root Children Pruner (3 (change to 2) levels, without repetition).
# During the loop phase,
# this strategy applies the Molecule Leaves Pruner, the High P-Value Branch Pruner (0.05), the Linear Branch Collapser Pruner,
# and the Zero Degree Vertex Pruner. No pruners are applied in the final phase post-loop.


def run_enrichment_analysis_plain_enrich_pruning_strategy(
    studyset_list: list[str],
    levels: int = 2,
    n: int = 0,
    p_value_threshold: float = 0.05,
    classification: str = "structural",
    print_results: bool = False,
    csv_output_path: str | None = None,
) -> tuple[dict, object]:
    """
    Run enrichment analysis with the Plain Enrichment Pruning Strategy.
    This strategy first applies the High Value Branch Pruner, the Linear Branch Collapser Pruner, and the Root Children Pruner in the pre-loop phase.
    Then, in the loop phase, it applies the High P-Value Branch Pruner, the Linear Branch Collapser Pruner, and the Zero Degree Vertex Pruner.

    Args:
        studyset_list: List of ChEBI class IDs (IRIs) to analyse.
        levels: Number of levels to prune from root (int, default 2). Only used for Root Children Pruner.
        n: Keep every n-th node in linear branches (int, default 0). Only used for Linear Branch Collapser Pruner.
        p_value_threshold: Threshold for high p-value pruning (float, default 0.05). Only used for High P-Value Branch Pruner.
        classification: Classification type for enrichment analysis ("structural", "functional", or "full").
        print_results: Print a results table to stdout before returning (bool, default False).
        csv_output_path: If given, write the enrichment results to this CSV path (default None).

    Returns:
        tuple: (results_dict, graph) where results_dict contains keys:
            "study_set" (item names), "removed_nodes" (pruned node names), and
            "enrichment_results" (node -> p-value mapping).
    """

    # Load common data files
    (
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        removed_leaves_csv,
        leaf_to_ancestors_map_file,
        parent_map_file,
    ) = _load_data_files()

    structural_leaf_ids = get_structural_leaf_ids(removed_leaves_csv)

    studyset_list = [normalize_id(cls) for cls in studyset_list]

    studyset_leaves = get_leaves(
        studyset_list,
        removed_leaves_csv,
        class_to_leaf_map,
        structural_leaf_ids,
    )
    print(f"Study set leaves: {describe(studyset_leaves)}")

    studyset_ancestors = get_ancestors_for_inputs(
        studyset_leaves,
        leaf_to_ancestors_map_file,
    )
    print(f"Study set ancestors: {describe(studyset_ancestors)}")
    print(f"Number of study set ancestors: {len(studyset_ancestors)}")

    all_removed_nodes = set()

    if classification in ["functional", "full"]:
        # Collect roles associated with leaves AND structural ancestors
        studyset_ancestors_roles = set()

        # Add roles from leaves
        for leaf in studyset_leaves:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(leaf, []))

        # Add roles from structural ancestors
        for class_to_check in studyset_ancestors:
            studyset_ancestors_roles.update(
                class_to_all_roles_map.get(class_to_check, []),
            )

        # Also include all ancestors of these roles (they will appear in the graph)
        with open(parent_map_file) as f:
            parent_map = json.load(f)

        roles_with_ancestors = set(studyset_ancestors_roles)
        to_process = list(studyset_ancestors_roles)

        while to_process:
            role = to_process.pop(0)
            for parent in parent_map.get(role, []):
                if parent not in roles_with_ancestors:
                    roles_with_ancestors.add(parent)
                    to_process.append(parent)

        studyset_ancestors_roles = roles_with_ancestors
        print(
            f"Number of roles (including all ancestors): {len(studyset_ancestors_roles)}",
        )
    else:
        studyset_ancestors_roles = set()

    enrichment_results = get_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles,
        structural_leaf_ids,
    )

    # print("Enrichment results:")
    # print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    # print("Enrichment results after Benjamini-Hochberg correction:")
    # print_enrichment_results(enrichment_results)

    pre_pruned_G = create_graph_with_roles_and_structures(
        studyset_leaves,
        studyset_ancestors,
        studyset_ancestors_roles,
        parent_map_file,
        class_to_all_roles_map,
        classification,
    )
    G = pre_pruned_G.copy()

    ## Pre-loop phase ##
    print("Starting pre-loop pruning phase.")
    G, removed_nodes = high_p_value_branch_pruner(
        G,
        enrichment_results,
        p_value_threshold,
    )
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by high p-value pruner: {describe(removed_nodes)}")

    G, removed_nodes = linear_branch_collapser_pruner_remove_less(G, n)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by linear branch pruner: {describe(removed_nodes)}")

    G, removed_nodes, _execution_count = root_children_pruner(
        G,
        levels,
        allow_re_execution=False,
        execution_count=0,
    )
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by root children pruner: {describe(removed_nodes)}")

    ## Loop phase ##
    print("Starting loop pruning phase.")
    # Count the number of nodes in G so it can be compared after each iteration
    size_before = G.number_of_nodes()
    size_after = size_before
    first_iteration = True
    iteration = 0

    # while the size changes, keep applying the loop phase pruners
    while size_after < size_before or first_iteration:
        size_before = size_after
        iteration += 1
        print(f"Loop iteration {iteration}")

        ## Recalculate corrected p-values

        # Remove pruned nodes from enrichment results
        current_enrichment = {
            cls: vals
            for cls, vals in enrichment_results.items()
            if cls not in all_removed_nodes and G.has_node(cls)
        }

        current_enrichment = benjamini_hochberg_fdr_correction(current_enrichment)

        G, removed_nodes = high_p_value_branch_pruner(
            G,
            current_enrichment,
            p_value_threshold,
        )
        all_removed_nodes.update(removed_nodes)
        print(f"Removed nodes by high p-value pruner: {describe(removed_nodes)}")

        G, removed_nodes = linear_branch_collapser_pruner_remove_less(G, n)
        all_removed_nodes.update(removed_nodes)
        print(f"Removed nodes by linear branch pruner: {describe(removed_nodes)}")

        G, removed_nodes = zero_degree_pruner(G)
        all_removed_nodes.update(removed_nodes)
        print(f"Removed nodes by zero-degree pruner: {describe(removed_nodes)}")

        size_after = G.number_of_nodes()
        first_iteration = False

    ## No final phase pruners ##

    final_enrichment = current_enrichment

    print(f"Number of removed nodes in total: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {
            id_to_name(cls): vals for cls, vals in final_enrichment.items()
        },
    }
    if print_results:
        print_enrichment_results(results["enrichment_results"])
    if csv_output_path:
        write_enrichment_csv(results["enrichment_results"], csv_output_path)

    return results, G


def run_enrichment_analysis_plain_enrich_pruning_strategy_from_smiles(
    smiles_list: list[str],
    use_parents: bool = False,
    **kwargs,
) -> tuple[dict, object, dict]:
    """Run run_enrichment_analysis_plain_enrich_pruning_strategy on SMILES strings.

    Each SMILES is resolved to ChEBI ID(s) via chebin.calculations.smiles_lookup
    before delegating to run_enrichment_analysis_plain_enrich_pruning_strategy; see
    that function for the remaining arguments (passed through as **kwargs) and
    return value.

    Args:
        smiles_list: List of SMILES strings to analyse.
        use_parents: If a SMILES can't be resolved to a direct ChEBI ID, fall back
            to a remote classification call and use its direct parent ChEBI IDs.

    Returns:
        tuple: (results_dict, graph, smiles_diagnostics) where smiles_diagnostics
            is {"unresolved_smiles": [...], "ambiguous_matches": [...]}.
    """
    studyset_list, unresolved_smiles, ambiguous_matches = smiles_list_to_studyset(
        smiles_list,
        use_parents=use_parents,
    )
    results, G = run_enrichment_analysis_plain_enrich_pruning_strategy(
        studyset_list,
        **kwargs,
    )
    return (
        results,
        G,
        {
            "unresolved_smiles": unresolved_smiles,
            "ambiguous_matches": ambiguous_matches,
        },
    )


if __name__ == "__main__":
    bonferroni_correct = False
    benjamini_hochberg_correct = True  # Used in Binche1

    root_children_prune = True
    levels = 2  # Number of levels to prune including root. 1 only prunes root, 2 prunes root, and it's direct neighbor, and so on.
    # allow_re_execution = False  # Currently not necessary. Whether the pruner can be executed multiple times on a given graph.
    # execution_count = 0  # Currently not necessary. Counter for the number of executions

    linear_branch_prune = True
    n = 2  # Keep only every n-th node in linear branches

    high_p_value_prune = True
    p_value_threshold = 0.05

    zero_degree_prune = True

    classification = "functional"  # "functional" or "structural" or "full"

    # For testing purposes, you can use the following study set of two compounds:

    studyset_list = [
        "http://purl.obolibrary.org/obo/CHEBI_77030",
        "http://purl.obolibrary.org/obo/CHEBI_79036",
    ]

    results = run_enrichment_analysis(
        studyset_list,
        bonferroni_correct=bonferroni_correct,
        benjamini_hochberg_correct=benjamini_hochberg_correct,
        root_children_prune=root_children_prune,
        levels=levels,
        linear_branch_prune=linear_branch_prune,
        n=n,
        high_p_value_prune=high_p_value_prune,
        p_value_threshold=p_value_threshold,
        zero_degree_prune=zero_degree_prune,
        classification=classification,
    )

    # run_enrichment_analysis_plain_enrich_pruning_strategy(studyset_list,
    #                         levels=2, # for root children pruner
    #                         n=2, # for linear branch pruner
    #                         p_value_threshold=0.05, # for high p-value pruner
    #                         classification="structural")
