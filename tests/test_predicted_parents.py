"""Tests for the rules on Chebifier-predicted parents (deepest only, leaf limit).

Runs against a toy hierarchy, so nothing here needs the data folder or network:

    root (CHEBI:1)
     └─ general (CHEBI:2, 200 leaves)
         ├─ specific (CHEBI:3, 10 leaves)
         └─ broad sibling (CHEBI:4, 500 leaves)
"""

import pytest

from chebin.calculations import predicted_parents, smiles_lookup
from chebin.calculations.predicted_parents import (
    LEAF_EXPANSION_LIMIT,
    select_expandable_parents,
)

IRI = "http://purl.obolibrary.org/obo/CHEBI_"
VALID_REMOTE_SMILES = "CCOC(=O)CCCCCCN1CCCC1"

pytestmark = pytest.mark.usefixtures("empty_local_maps")


@pytest.fixture(autouse=True)
def toy_hierarchy(monkeypatch):
    parent_map = {
        f"{IRI}2": [f"{IRI}1"],
        f"{IRI}3": [f"{IRI}2"],
        f"{IRI}4": [f"{IRI}2"],
    }
    leaf_counts = {f"{IRI}1": 10_000, f"{IRI}2": 200, f"{IRI}3": 10, f"{IRI}4": 500}
    monkeypatch.setattr(
        predicted_parents,
        "_get_parent_map_and_leaf_counts",
        lambda: (parent_map, leaf_counts),
    )
    monkeypatch.setattr(predicted_parents, "_depth_memo", {})


def test_limit_matches_the_background_builder():
    assert LEAF_EXPANSION_LIMIT == 150


def test_only_the_deepest_parent_is_kept():
    assert select_expandable_parents(["CHEBI:2", "CHEBI:3"]) == (["CHEBI:3"], [])


def test_a_tied_deepest_parent_over_the_limit_is_skipped():
    assert select_expandable_parents(["CHEBI:3", "CHEBI:4"]) == (
        ["CHEBI:3"],
        ["CHEBI:4"],
    )


def test_a_lone_parent_over_the_limit_is_skipped():
    assert select_expandable_parents(["CHEBI:2"]) == ([], ["CHEBI:2"])


def test_no_parents():
    assert select_expandable_parents([]) == ([], [])


def _stub_chebifier(monkeypatch, parent_ids):
    def post(url, payload, label):
        if url == smiles_lookup.CHEBIFIER_DETAILS_URL:
            return {
                "models": {"ChEBI Lookup": {"chebi_ids": [], "highlights": []}},
            }, None
        return {"direct_parents": [[[pid, "name"] for pid in parent_ids]]}, None

    monkeypatch.setattr(smiles_lookup, "_post_chebifier_json", post)


def test_lookup_uses_only_the_expandable_parents(monkeypatch):
    _stub_chebifier(monkeypatch, ["2", "3", "4"])
    chebi_ids, was_resolved, _, failure_reason, too_general = (
        smiles_lookup.convert_smiles_to_chebi(VALID_REMOTE_SMILES, use_parents=True)
    )
    assert chebi_ids == ["CHEBI:3"]
    assert was_resolved is True
    assert failure_reason is None
    assert too_general == ["CHEBI:4"]


def test_lookup_reports_a_structure_whose_parents_are_all_too_general(monkeypatch):
    """Not "no ChEBI match": parents were found, they just couldn't be used."""
    _stub_chebifier(monkeypatch, ["4"])
    chebi_ids, was_resolved, _, failure_reason, too_general = (
        smiles_lookup.convert_smiles_to_chebi(VALID_REMOTE_SMILES, use_parents=True)
    )
    assert chebi_ids == []
    assert was_resolved is False
    assert failure_reason == "parents_too_general"
    assert too_general == ["CHEBI:4"]
