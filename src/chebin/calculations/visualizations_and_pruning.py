import json
import os
import xml.etree.ElementTree as ET
from collections import Counter
from math import inf

import networkx as nx

from chebin.calculations.chebi_ids import to_chebi_curie
from chebin.calculations.data_files import ID_TO_NAME_MAP
from chebin.config import require_data_path


def _id_to_name_map_file() -> str:
    """Resolved per call, so set_data_dir() still applies after import."""
    return require_data_path(ID_TO_NAME_MAP)


# (mtime, map) of the last loaded name map. id_to_name is called once per graph
# node, per tested class and per result entry -- thousands of times per analysis --
# so re-reading this 17 MB / ~225k-entry file per call dominated runtime. Cached on
# mtime rather than outright so the monthly data/ refresh (finalize_folder_structure)
# is picked up without restarting the server.
_id_to_name_cache: tuple[float, dict] | None = None


def _load_id_to_name_map(path: str | None = None) -> dict:
    """Return the ChEBI id->name map, re-reading it only when the file changes."""
    global _id_to_name_cache

    if path is None:
        path = _id_to_name_map_file()
    mtime = os.path.getmtime(path)
    if _id_to_name_cache is None or _id_to_name_cache[0] != mtime:
        with open(path) as f:
            _id_to_name_cache = (mtime, json.load(f))
    return _id_to_name_cache[1]


def id_to_name(class_id: str) -> str:
    """
    Convert ChEBI class IRI to human-readable name with ID.

    Uses the cached name mapping and returns formatted string
    "Name (CHEBI:ID)" or just the ID if name not found.

    Args:
        class_id: Class IRI (e.g., http://purl.obolibrary.org/obo/CHEBI_12345).

    Returns:
        str: Formatted name like "ascorbic acid (CHEBI:15377)" or just class_id.
    """
    id_to_name_map = _load_id_to_name_map()

    prefix = "http://purl.obolibrary.org/obo/"
    if class_id.startswith(prefix):
        # remove prefix
        class_id = class_id.replace(prefix, "")
    name = id_to_name_map.get(class_id)
    display_id = class_id.replace("_", ":", 1)
    return f"{name} ({display_id})" if name else display_id


def strip_prefix(class_id: str) -> str:
    """
    Remove OBO namespace prefix from class IRI.

    Converts http://purl.obolibrary.org/obo/CHEBI_12345 to CHEBI_12345.

    Args:
        class_id: Class IRI (with or without prefix).

    Returns:
        str: Class ID without the OBO prefix.
    """
    prefix = "http://purl.obolibrary.org/obo/"
    if class_id.startswith(prefix):
        return class_id.replace(prefix, "")
    return class_id


#####################################
# Forming graph
#####################################

# def find_paths_to_root_old(ontology, start_class):
#     paths = []

#     def dfs(current_class, current_path):
#         superclasses = ontology.get_superclasses(current_class)
#         superclasses = [s for s in superclasses if s not in current_path] # Remove circular references

#         if not superclasses: # Reached root
#             paths.append(current_path)
#             return

#         for superclass in superclasses:
#             dfs(superclass, current_path + [superclass]) # Appends superclass to path

#     dfs(start_class, [start_class])
#     return paths


def find_paths_to_root_with_map(
    start_class: str,
    parents_map: dict[str, list[str]],
) -> list[list[str]]:
    """
    Find all paths from a class to root(s) in the ontology hierarchy.

    Uses depth-first search to explore parent relationships until reaching
    a root class (one with no parents in the map).

    Args:
        start_class: Starting class IRI.
        parents_map: Dictionary mapping class IRIs to lists of parent IRIs.

    Returns:
        list: List of paths, where each path is a list of class IRIs from
            start_class to a root class (inclusive).
    """
    paths: list[list[str]] = []

    def dfs(current_class: str, current_path: list[str]) -> None:
        # if the class has no parents in the map, it's a root
        parents = parents_map.get(current_class, [])
        if not parents:
            paths.append(current_path)
            return

        for parent in parents:
            dfs(parent, current_path + [parent])

    dfs(start_class, [start_class])
    return paths


# Not used anymore
def find_paths_to_root_with_ontology(
    ontology,
    start_class: str,
    leaf_to_parents_json_file: str = "data/removed_leaf_classes_to_direct_parents_map.json",
) -> list[list[str]]:
    """Find all paths from start_class to ontology roots."""
    paths: list[list[str]] = []

    # Unified DFS (works for both ontology and leaf parents)
    def dfs(current_class: str, current_path: list[str]) -> None:
        superclasses = ontology.get_superclasses(current_class)
        superclasses = [s for s in superclasses if s not in current_path]

        if not superclasses:  # reached root
            paths.append(current_path)
            return

        for superclass in superclasses:
            dfs(superclass, current_path + [superclass])

    # -----------------------------------------
    # CASE 1: class is a leaf → load parent map
    # -----------------------------------------
    with open(leaf_to_parents_json_file) as f:
        leaf_to_parents = json.load(f)

    if start_class in leaf_to_parents:
        direct_parents = leaf_to_parents[start_class]
        if not isinstance(direct_parents, list):
            direct_parents = [direct_parents]

        # Start DFS from each parent, but include the leaf in the initial path
        for parent in direct_parents:
            dfs(parent, [start_class, parent])

        return paths

    else:
        # -----------------------------------------
        # CASE 2: class exists in filtered ontology (is not a leaf)
        # -----------------------------------------

        dfs(start_class, [start_class])
        return paths


def get_name(chebi_ontology: str, iri: str) -> str | None:
    """Get human-readable name for a class IRI from OWL ontology."""
    tree = ET.parse(chebi_ontology)
    root = tree.getroot()
    ns = {
        "owl": "http://www.w3.org/2002/07/owl#",
        "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    }

    # Find all OWL classes
    for cls in root.findall("owl:Class", ns):
        about = cls.attrib.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about")
        if about == iri:
            # print(f"Found class for IRI: {iri}")
            # Find rdfs:label element
            label_elem = cls.find("rdfs:label", ns)
            if label_elem is not None and label_elem.text is not None:
                return label_elem.text.strip()

    return None  # if not found


def create_graph_from_paths(paths: list) -> nx.DiGraph:
    G = nx.DiGraph()
    for path in paths:
        for i in range(len(path) - 1):
            G.add_edge(path[i], path[i + 1])

    return G


# # Not used anymore. If to be used, potentially remove color_map parameter
# def create_graph_from_ontology(classes, classification, color_map = ['#FFB6C1', "#F44280", "#AA83A7", "#83163A", "#E63FE6", '#FFA07A', '#FF69B4'], max_n_leaf_classes=inf):
#     G = nx.DiGraph()
#     j = 0

#     if classification == "structural":
#         ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl")
#     elif classification == "functional":
#         ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl")
#     else:
#         ontology = load_chebi()

#     for i, cls in enumerate(classes):
#         print(f"Adding to graph... Starting node: {cls}")
#         paths = find_paths_to_root(ontology, cls)
#         H = create_graph_from_paths(paths)
#         # Color the nodes of H
#         color = color_map[j % len(color_map)]
#         nx.set_node_attributes(H, color, 'color')

#         G = nx.compose(G, H)  # Combine graphs
#         j += 1
#         if j >= max_n_leaf_classes:
#             break

#         print(f"Total number of starting leaf classes processed in graph: {j}")
#     return G


def create_graph_from_map(
    classes: list[str],
    parent_map_json_file: str,
    max_n_leaf_classes: float = inf,
) -> nx.DiGraph:
    """
    Create a directed graph of class hierarchy from a parent map JSON file.

    Constructs a networkx DiGraph where nodes are classes and edges represent
    parent-child relationships. Uses paths from each class to root(s).

    Args:
        classes: Iterable of starting class IRIs (typically leaves).
        parent_map_json_file: Path to JSON file with parent mapping.
        max_n_leaf_classes: Maximum number of starting classes to process
            (for memory efficiency on large datasets).

    Returns:
        networkx.DiGraph: Directed graph of the class hierarchy.
    """

    with open(parent_map_json_file) as f:
        parents_map = json.load(f)

    G: nx.DiGraph = nx.DiGraph()
    j = 0

    for cls in classes:
        if j % 10 == 0:
            print(f"Processing class {j + 1}/{len(classes)}")

        paths = find_paths_to_root_with_map(cls, parents_map)

        for path in paths:
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]

                if not G.has_edge(u, v):
                    G.add_edge(u, v)

                    # add label once, when node first appears
                    if "label" not in G.nodes[u]:
                        G.nodes[u]["label"] = id_to_name(u)
                    if "label" not in G.nodes[v]:
                        G.nodes[v]["label"] = id_to_name(v)

        j += 1
        if j >= max_n_leaf_classes:
            break

    print(f"Total number of starting leaf classes processed in graph: {j}")
    return G


def create_graph_with_roles_and_structures(
    studyset_leaves: list[str],
    structural_ancestors: list[str] | set[str],
    enriched_roles: list[str] | set[str],
    parent_map_file: str,
    class_to_all_roles_map: dict,
    classification: str,
) -> nx.DiGraph:

    if classification == "structural" or classification == "full":
        # Build structural graph
        G = create_graph_from_map(studyset_leaves, parent_map_file)

    elif classification == "functional":
        # Build role hierarchy only
        G = nx.DiGraph()

    if classification == "functional" or classification == "full":
        # Build role hierarchy
        with open(parent_map_file) as f:
            parents_map = json.load(f)

        for role in enriched_roles:
            # Build path from role to root
            paths = find_paths_to_root_with_map(role, parents_map)
            for path in paths:
                for i in range(len(path) - 1):
                    u, v = path[i], path[i + 1]

                    if not G.has_edge(u, v):
                        G.add_edge(u, v)

                        # add label once, when node first appears
                        if "label" not in G.nodes[u]:
                            G.nodes[u]["label"] = id_to_name(u)
                        if "label" not in G.nodes[v]:
                            G.nodes[v]["label"] = id_to_name(v)

    return G


# # Similar to the above but doesn't first create separate graphs gor each class
# def create_graph_from_map_original(classes, parent_map_json_file, max_n_leaf_classes=inf):

#     with open(parent_map_json_file, "r") as f:
#         parents_map = json.load(f)

#     G = nx.DiGraph()
#     j = 0


#     for cls in classes:
#         if j % 10 == 0:
#             print(f"Proceeesing class {j+1}/{len(classes)}")

#         paths = find_paths_to_root_with_map(cls, parents_map)
#         H = create_graph_from_paths(paths)

#         # Add labels for Cytospace compatibility
#         label_dict = {node: id_to_name(node) for node in H.nodes()}
#         nx.set_node_attributes(H, label_dict, 'label')

#         G = nx.compose(G, H)  # Combine graphs

#         j += 1
#         if j >= max_n_leaf_classes:
#             break

#         print(f"Total number of starting leaf classes processed in graph: {j}")
#     return G


# def draw_graph(G, graphing_layout, title):
#     if graphing_layout == "default":
#         pos = None  # Default layout
#     elif graphing_layout == "kamada_kawai":
#         pos = nx.kamada_kawai_layout(G)
#     elif graphing_layout == "spectral":
#         pos = nx.spectral_layout(G)
#     elif graphing_layout == "layer_based":
#         # Calculate depth of each node from root
#         roots = [n for n, d in G.in_degree() if d == 0]
#         if roots:
#             # Assign layer based on shortest path from root
#             layers = {}
#             for node in G.nodes():
#                 min_dist = float('inf')
#                 for root in roots:
#                     if nx.has_path(G, root, node):
#                         dist = nx.shortest_path_length(G, root, node)
#                         min_dist = min(min_dist, dist)
#                 layers[node] = min_dist if min_dist != float('inf') else 0

#             # Set subset attribute for multipartite layout
#             nx.set_node_attributes(G, layers, 'subset')
#             pos = nx.multipartite_layout(G, subset_key='subset', align='horizontal')
#     else:
#         print(f"Unknown graphing layout: {graphing_layout}. Using default.")
#         pos = None  # Default layout


#     plt.figure(figsize=(20, 10))

#     # Draw nodes with their assigned colors
#     node_colors = [G.nodes[n].get("color") for n in G.nodes()]

#     nx.draw(G, pos, with_labels=True, node_size=500, node_shape='s', font_size=8, font_weight='bold', node_color=node_colors, arrows=True,
#             arrowsize=12, edge_color='black', alpha=1)

#     plt.title(title, fontsize=12)
#     plt.show()

#####################################
# Pruning strategies
#####################################

# RootChildrenPruner:
# This pruner deletes roots and their children up to a defined level.
### 2 levels remove two levels: root and its direct children.
### Different from Binche1 which only removes root and 2 of its direct children.


def delete_children(node, G, next_level, removed_nodes):

    # Traverse down to the specified level and remove nodes
    if next_level > 0:
        next_level -= 1
        children = list(
            G.predecessors(node),
        )  # In DiGraph, predecessors are children in this case
        for child in children:
            delete_children(child, G, next_level, removed_nodes)

        # Record note removal
        G.remove_node(node)
        removed_nodes.add(node)


def root_children_pruner(G, levels, allow_re_execution=False, execution_count=0):
    """
    Remove nodes from levels 1 to `levels` of the root (not including root itself).

    Args:
        G (networkx.DiGraph): Ontology graph
        levels (int): Number of levels to prune from root
        allow_re_execution (bool): If False, only executes on first call (count=0)
                                   If True, can execute multiple times
        execution_count (int): Counter tracking how many times this pruner has run

    Returns:
        tuple: (G, removed_nodes, execution_count) where removed_nodes is the
               set of removed node IDs, and execution_count is incremented if
               the pruner executed.

    Note:
        The allow_re_execution and execution_count parameters track whether
        pruners can run multiple times in a pipeline. Currently these are used
        but the logic is basic; potential future optimization to streamline.
    """
    removed_nodes = set()
    if allow_re_execution or execution_count == 0:
        roots = [n for n, d in G.out_degree() if d == 0]
        for root in roots:
            delete_children(root, G, levels, removed_nodes)
        execution_count += 1
    return G, removed_nodes, execution_count


# Linear branch collapser pruner - remove fewer nodes:
# This pruner collapses linear branches in the ontology graph by removing
# fewer intermediate nodes in branches where each node has exactly one child,
# effectively connecting every n:th node in such branches directly.
# Keeps every n:th node in linear branches.


# new version
def process_branch_remove_less(head, node, G, n, removed_nodes):
    # Check if node still exists (might have been removed in another branch)
    if not G.has_node(node):
        return G

    branch_nodes = []
    current_node = node
    children = list(G.predecessors(current_node))

    while len(children) == 1:
        branch_nodes.append(current_node)
        current_node = children[0]
        # Check if current_node still exists before getting its predecessors
        if not G.has_node(current_node):
            break
        children = list(G.predecessors(current_node))

    # print(f"Current node: {current_node}, Children: {children}")
    last_node = current_node

    # Capture children BEFORE modifying the graph
    if G.has_node(last_node):
        children_before = list(G.predecessors(last_node))
    else:
        children_before = []

    if len(branch_nodes):
        # Determine nodes to keep
        if n == 0:
            # Remove all intermediate nodes
            nodes_to_keep = [head, last_node]
        else:
            # Keep every n:th node in the branch
            nodes_to_keep = [head]
            for index, branch_node in enumerate(branch_nodes, start=1):
                if index % n == 0:
                    nodes_to_keep.append(branch_node)
            nodes_to_keep.append(last_node)

            # print(f"Nodes to keep in branch: {nodes_to_keep}")

        for branch_node in branch_nodes:
            if branch_node not in nodes_to_keep:
                removed_nodes.add(branch_node)
                G.remove_node(branch_node)

        for i in range(len(nodes_to_keep) - 1):
            if (
                G.has_node(nodes_to_keep[i])
                and G.has_node(nodes_to_keep[i + 1])
                and not G.has_edge(nodes_to_keep[i + 1], nodes_to_keep[i])
            ):
                G.add_edge(nodes_to_keep[i + 1], nodes_to_keep[i])

    # Recurse only over children that still exist in the graph
    for child in children_before:
        if G.has_node(child):  # Check before recursing
            process_branch_remove_less(last_node, child, G, n, removed_nodes)

    return G


def linear_branch_collapser_pruner_remove_less(G, n):
    removed_nodes = set()
    roots = [n for n, d in G.out_degree() if d == 0]
    for root in roots:
        # print(f"Processing root: {root}")
        direct_children = list(G.predecessors(root))
        # print(f"Direct children of root {root}: {direct_children}")
        for child in direct_children:
            process_branch_remove_less(root, child, G, n, removed_nodes)
    return G, removed_nodes


# High P-Value Branch Pruner:
# Removes branches from the graph components that contain only vertices with a
# p-value greater than 0.05. If the branch inspected has at least one node
# with a p-value below the threshold, the branch is kept.


def _new_p_value_pruner_stats():
    """Accumulator for conditions worth reporting once per run rather than per visit."""
    return {
        "no_pvalue_visits": 0,
        "no_pvalue_leaf": set(),
        "no_pvalue_non_leaf": set(),
        "already_removed_visits": 0,
        "already_removed_nodes": set(),
    }


def _log_p_value_pruner_stats(stats):
    """Summarise the pruner walk in a few lines.

    These conditions used to print once per node *visit*. Because the walk below
    is not memoised and ChEBI is a DAG with heavy multiple inheritance, a node is
    re-entered once per distinct root path, which produced ~400k near-identical
    lines in a single real run. Aggregating mirrors graph_to_cytospace_json.
    """
    leaves = stats["no_pvalue_leaf"]
    non_leaves = stats["no_pvalue_non_leaf"]

    if leaves or non_leaves:
        print(
            f"High p-value pruner: {len(leaves) + len(non_leaves)} nodes had no p-value "
            f"across {stats['no_pvalue_visits']} visits; {len(leaves)} are untested "
            f"study-set leaves (expected, they are never assigned p-values).",
        )

    # A node with children should have been tested, so this one is worth seeing.
    if non_leaves:
        sample = sorted(non_leaves)[:10]
        print(
            f"High p-value pruner: WARNING {len(non_leaves)} non-leaf node(s) had no "
            f"p-value, which is unexpected. First {len(sample)}: {sample}",
        )

    if stats["already_removed_nodes"]:
        print(
            f"High p-value pruner: revisited {len(stats['already_removed_nodes'])} "
            f"already-removed node(s) across {stats['already_removed_visits']} visits.",
        )


def high_p_value_branch_pruner(G, p_value_dict, p_value_threshold=0.05):
    removed_nodes = set()
    roots = [n for n, d in G.out_degree() if d == 0]
    stats = _new_p_value_pruner_stats()

    for root in roots:
        size_before = G.number_of_nodes()
        process_node_for_p_value_pruner(
            root,
            G,
            p_value_dict,
            p_value_threshold,
            removed_nodes,
            stats,
        )
        size_after = G.number_of_nodes()
        # Keep pruning until nothing more is removed
        while size_after < size_before:
            size_before = size_after
            process_node_for_p_value_pruner(
                root,
                G,
                p_value_dict,
                p_value_threshold,
                removed_nodes,
                stats,
            )
            size_after = G.number_of_nodes()

    _log_p_value_pruner_stats(stats)

    return G, removed_nodes


def process_node_for_p_value_pruner(
    node,
    G,
    p_value_dict,
    p_value_threshold,
    removed_nodes,
    stats=None,
):  # Returns a boolean. Tells whether node or any descendant has p-value below threshold

    # Check if node still exists (might have been removed in another branch)
    if not G.has_node(node):
        if stats is not None:
            stats["already_removed_visits"] += 1
            stats["already_removed_nodes"].add(node)
        return False

    has_good_descendant = False  # whether any descendant has p-value below threshold
    children = list(G.predecessors(node))
    nodes_to_remove = []
    for child in children:
        if process_node_for_p_value_pruner(
            child,
            G,
            p_value_dict,
            p_value_threshold,
            removed_nodes,
            stats,
        ):
            has_good_descendant = True
        else:
            nodes_to_remove.append(child)

    # Node's own p-value check
    # Get correct p-value for the node
    if "p_value_corrected" in p_value_dict.get(node, {}):
        node_p_value = p_value_dict.get(node, {}).get("p_value_corrected")
    else:
        node_p_value = p_value_dict.get(node, {}).get("p_value", None)

    if node_p_value is None:
        # Expected for study-set leaves: the graph is built from leaf->root paths,
        # but only ancestors are tested, so leaves never receive a p-value. Counted
        # instead of printed -- see _log_p_value_pruner_stats. A node that still has
        # children should have been tested, so it is tracked separately as unexpected.
        if stats is not None:
            stats["no_pvalue_visits"] += 1
            if children:
                stats["no_pvalue_non_leaf"].add(node)
            else:
                stats["no_pvalue_leaf"].add(node)
        # Report False: this node contributes no evidence of significance, which is
        # what the return value means to the caller. It previously returned True to
        # express "don't delete me", but the caller reads that as "this branch is
        # significant", so a single untested leaf shielded its entire ancestor chain
        # from pruning. Returning early (before the removal branch below) already
        # keeps the node itself; any leaf left stranded by its pruned parent is
        # cleaned up by zero_degree_pruner.
        return False

    if node_p_value <= p_value_threshold:
        has_good_descendant = True
    elif children == [] or (len(nodes_to_remove) == len(children)):
        # No good descendants and node itself is bad → mark for removal
        removed_nodes.add(node)
        if G.has_node(node):  # Check before removing
            G.remove_node(node)
        return False  # Whole branch is bad

    return has_good_descendant


def zero_degree_pruner(G: nx.DiGraph) -> tuple[nx.DiGraph, list]:
    to_remove = []
    removed = []

    for node in G.nodes():
        if G.degree(node) == 0:
            to_remove.append(node)

    for node in to_remove:
        G.remove_node(node)
        removed.append(node)

    return G, removed


#####################################
# Converting NetworkX graph to Cytoscape compatible format
#####################################
def clean_label(label: str) -> str:
    """Changes label from 'Name (CHEBI_ID)' to 'Name'"""
    return label.split(" (")[0] if " (" in label else label


def extract_chebi_id(label: str) -> str | None:
    """Finds CHEBI_XXXX inside parentheses"""
    if " (" in label and "CHEBI_" in label:
        return label.split(" (")[1][:-1]  # Extract CHEBI_ID without parentheses
    return label


def graph_to_cytospace_json(
    G,
    output_file,
    enrichment_results=None,
    include_untested_leaves=False,
):
    """Write the graph as Cytoscape JSON to `output_file`.

    Thin wrapper around :func:`graph_to_cytoscape_dict`; see there for the meaning of
    `include_untested_leaves`.
    """
    data = graph_to_cytoscape_dict(
        G,
        enrichment_results=enrichment_results,
        include_untested_leaves=include_untested_leaves,
    )

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as f:
        json.dump(data, f, indent=4)


def graph_to_cytoscape_dict(
    G,
    enrichment_results=None,
    include_untested_leaves=False,
):
    """Build the Cytoscape elements structure for the graph.

    include_untested_leaves: study-set leaf classes are graph nodes (the graph is
    built from leaf->root paths) but are never tested, so they can never carry a
    p-value and are not part of the significance view. They typically make up the
    large majority of nodes -- ~86% in a real run -- so shipping them forces the
    browser to lay out thousands of nodes it will never colour. They are excluded
    by default; pass True to keep them for debugging. This only affects what is
    drawn: pruning and all p-values are already final by this point, and the full
    study set is still reported in enrichment_results["study_set"].

    Returns:
        dict: ``{"elements": [...]}`` ready for Cytoscape.
    """
    data = {"elements": []}

    # Extract the nested enrichment results
    enr_dict = (
        enrichment_results.get("enrichment_results") if enrichment_results else None
    )
    study_set_labels = (
        set(enrichment_results.get("study_set", [])) if enrichment_results else set()
    )
    missing_pvalue_nodes = []
    skipped_leaves = set()

    # Nodes
    for node, attrs in G.nodes(data=True):
        label = attrs.get("label", node)
        short_label = clean_label(label)
        node_data = {
            "id": node,
            "label": label,
            "short_label": short_label,
            "color": attrs.get("color", "#706C6C"),
            # The node's identity as ChEBI itself writes it, so the graph view can
            # link straight to the entry page without re-parsing the IRI or scraping
            # the id back out of the label. None for the rare node that is not a
            # ChEBI class; the view leaves those unlinked.
            "chebi_id": to_chebi_curie(node),
        }
        # Add enrichment results to nodes if provided

        if enr_dict is not None and label in enr_dict:
            enr = enr_dict[label]
            node_data["p_value"] = enr.get("p_value")
            node_data["p_value_corrected"] = enr.get("p_value_corrected")
            node_data["p_value_reason"] = "present"
        else:
            if enr_dict is None:
                reason = "no_enrichment_results_dict"
            elif label in study_set_labels:
                reason = "study_set_leaf_not_tested"
            else:
                reason = "node_not_in_enrichment_results"
            node_data["p_value_reason"] = reason
            missing_pvalue_nodes.append((label, reason))

            # Untested study-set leaves are not part of the significance view.
            # Note this deliberately keeps "node_not_in_enrichment_results" nodes,
            # which are unexpected and worth seeing.
            if not include_untested_leaves and reason == "study_set_leaf_not_tested":
                skipped_leaves.add(node)
                continue

        data["elements"].append({"data": node_data})

    # enrichment_results entries contain class metadata, counts, odds ratio,
    # and p-values; the JSON structure is built above.

    # Edges (skip any edge touching an omitted leaf, so no edge dangles)
    for source, target in G.edges():
        if source in skipped_leaves or target in skipped_leaves:
            continue
        data["elements"].append(
            {
                "data": {
                    "id": f"{source}_to_{target}",
                    "source": source,
                    "target": target,
                },
            },
        )

    if skipped_leaves:
        print(
            f"Graph view: omitted {len(skipped_leaves)} untested study-set leaf "
            f"node(s) from the graph JSON (they never receive p-values); "
            f"{G.number_of_nodes() - len(skipped_leaves)} node(s) written.",
        )

    if missing_pvalue_nodes:
        total_nodes = G.number_of_nodes()
        missing_count = len(missing_pvalue_nodes)
        reason_counts = Counter(reason for _, reason in missing_pvalue_nodes)
        print(
            f"Graph p-value coverage: {total_nodes - missing_count}/{total_nodes} nodes with p-values; "
            f"{missing_count} nodes with N/A",
        )
        print(f"Graph N/A reasons: {dict(reason_counts)}")
        for label, reason in missing_pvalue_nodes[:15]:
            print(f"Graph N/A node: {label} | reason={reason}")

    return data


# Example usage:
# start_class = "http://purl.obolibrary.org/obo/CHEBI_33675"
# start_time = time.time()
# ontology = load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl")
# paths = find_paths_to_root(ontology, start_class)
# end_time = time.time()
# print(f"Time taken using ontology: {end_time - start_time} seconds")
# print("Using ontology:")
# for path in paths:
#     print(" -> ".join(path))
# start_time = time.time()
# paths = find_paths_to_root_with_map("data/chebi_parent_map.json", start_class)
# end_time = time.time()
# print(f"Time taken using parent map: {end_time - start_time} seconds")
# print("Using parent map:")
# for path in paths:
#     print(" -> ".join(path))
#
# ---- Usage ----
# Variables
# levels = 2 # Number of levels to prune from root. 1 only prunes root, and it's direct neighbor, and so on.
# allow_re_execution = False  # True or False. whether the pruner can be executed multiple times on a given graph.
# execution_count = 0  # Counter for the number of executions
#
# OBS: first run cell to create G, then run the pruner function below.
# pruned_G = G.copy()
# pruned_G, execution_count = root_children_pruner(pruned_G, levels, allow_re_execution, execution_count)
#
# graphing_layout = "kamada_kawai" # options: "default", "kamada_kawai", "spectral", "layer_based"
# draw_graph(pruned_G, graphing_layout, f"Ontology graph pruned {levels} levels from root(s)")
