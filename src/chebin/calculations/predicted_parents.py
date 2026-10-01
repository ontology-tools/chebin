"""Rules for which Chebifier-predicted parent classes may stand in for a structure.

When a structure has no ChEBI entry of its own, Chebifier predicts parent classes
for it, which are then expanded to all their leaf descendants. Left unchecked, a
generic prediction (e.g. "organooxygen compound", ~115k leaves) swamps whatever
it is added to. The same two rules therefore apply wherever predicted parents are
used -- building the restricted backgrounds and resolving study-set SMILES on the
website -- so a structure ends up represented the same way on both sides:

1. Of several predicted parents, keep only the deepest in the hierarchy.
2. Don't expand a class with more than :data:`LEAF_EXPANSION_LIMIT` leaves.
"""

from __future__ import annotations

import json

from chebin.calculations.data_files import CLASS_TO_LEAF_MAP, PARENT_MAP
from chebin.config import require_data_path

#: A class with more leaf descendants than this is not expanded, to keep
#: high-level classes from inflating a background or study set.
LEAF_EXPANSION_LIMIT = 150

_CHEBI_IRI_PREFIX = "http://purl.obolibrary.org/obo/CHEBI_"


def _calculate_depth_to_root(chebi_iri, parent_map, memo=None):
    """Recursively calculate depth (path length) from a node to the root.

    Uses memoization to avoid recalculating already-seen nodes.
    If multiple parents exist, returns the maximum depth among them.
    """
    if memo is None:
        memo = {}
    if chebi_iri in memo:
        return memo[chebi_iri]

    parents = parent_map.get(chebi_iri, [])
    if not parents:
        memo[chebi_iri] = 0
        return 0

    # depth = 1 + maximum depth among all parents
    depth = 1 + max(_calculate_depth_to_root(p, parent_map, memo) for p in parents)
    memo[chebi_iri] = depth
    return depth


def filter_chebifier_parents(parent_chebis, chebi_parent_map, memo=None):
    """If Chebifier finds several parent classes for a given class (that does not have its own CHEBI ID),
    we want to only keep the parent(s) furthest down the hierarchy, i.e. the one(s) with the longest path to the root.
    This is to avoid inflating the background with very high-level classes.

    Input: list of parent CHEBI IDs (strings), and the CHEBI parent-child map --
           either a ready dict or a path to its JSON file. Callers iterating over
           many entities should load the JSON once and pass the dict, since this
           function is typically called once per multi-ID entity.
    Output: list of parent CHEBI IDs (strings) that are furthest down the hierarchy.
            should only be one parent unless there is a tie
    """

    if isinstance(chebi_parent_map, str):
        with open(chebi_parent_map, encoding="utf-8") as f:
            chebi_parent_map = json.load(f)

    # calculate path length to root for each parent CHEBI ID using memoization
    if memo is None:
        memo = {}
    parent_path_lengths = {
        p: _calculate_depth_to_root(p, chebi_parent_map, memo) for p in parent_chebis
    }

    # find the maximum path length
    max_path_length = max(parent_path_lengths.values())

    # choose the parent(s) with the maximum path length
    chosen_parents = [
        parent
        for parent, path_length in parent_path_lengths.items()
        if path_length == max_path_length
    ]

    return chosen_parents


_parent_map = None
_leaf_counts = None
_depth_memo: dict[str, int] = {}


def _get_parent_map_and_leaf_counts():
    """Lazily load and memoize the parent map and per-class leaf counts.

    Only the counts are kept from the (large) class-to-leaf map, since the limit
    check needs nothing else. Loaded on first use so that importing this module
    doesn't require the data folder to be configured yet.
    """
    global _parent_map, _leaf_counts
    if _parent_map is None or _leaf_counts is None:
        with open(require_data_path(PARENT_MAP), encoding="utf-8") as f:
            _parent_map = json.load(f)
        with open(require_data_path(CLASS_TO_LEAF_MAP), encoding="utf-8") as f:
            _leaf_counts = {cls: len(leaves) for cls, leaves in json.load(f).items()}
    return _parent_map, _leaf_counts


def select_expandable_parents(parent_ids):
    """Apply both rules to one structure's predicted parents.

    parent_ids are ``CHEBI:123`` CURIEs as Chebifier's results are turned into.
    Returns (kept, skipped_too_general), both as CURIEs: kept are the deepest
    parents within the leaf limit; skipped_too_general are the deepest parents
    left out for exceeding it. Parents dropped for not being the deepest are in
    neither list -- the background builder doesn't report those either.
    """
    if not parent_ids:
        return [], []
    parent_map, leaf_counts = _get_parent_map_and_leaf_counts()

    iris = [_CHEBI_IRI_PREFIX + cid.split(":", 1)[1] for cid in parent_ids]
    deepest = (
        filter_chebifier_parents(iris, parent_map, _depth_memo)
        if len(iris) > 1
        else iris
    )

    kept, skipped = [], []
    for iri in deepest:
        curie = "CHEBI:" + iri.removeprefix(_CHEBI_IRI_PREFIX)
        if leaf_counts.get(iri, 0) > LEAF_EXPANSION_LIMIT:
            skipped.append(curie)
        else:
            kept.append(curie)
    return kept, skipped
