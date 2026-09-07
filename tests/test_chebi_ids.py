"""Tests for chebin.calculations.chebi_ids -- the ChEBI ID formats users type.

The website only ever accepted ``CHEBI:17079``; anything else a user pasted
(``chebi:17079``, ``CHEBI ID: 17079``, a bare ``17079``) reached the enrichment
analysis as a node matching nothing, so the run silently came back empty. These
cover the forms that now have to be folded together, including the ones the
website's line splitter used to tear apart.
"""

import pytest

from chebin.calculations.chebi_ids import (
    CHEBI_IRI_PREFIX,
    collapse_chebi_prefixes,
    looks_like_chebi_id,
    parse_chebi_number,
    to_chebi_curie,
    to_chebi_iri,
)
from chebin.calculations.fishers_calculations import normalize_id
from chebin.calculations.smiles_lookup import is_smiles

ACCEPTED_FORMS = [
    "CHEBI:17079",
    "chebi:17079",
    "ChEBI:17079",
    "CHEBI_17079",
    "chebi_17079",
    "CHEBI 17079",
    "CHEBI-17079",
    "CHEBI17079",
    "CHEBI ID: 17079",
    "CHEBI ID:17079",
    "chebi id 17079",
    "17079",
    "  CHEBI:17079  ",
    '"CHEBI:17079"',
    "CHEBI:0017079",
    "http://purl.obolibrary.org/obo/CHEBI_17079",
    "https://identifiers.org/CHEBI:17079",
]


@pytest.mark.parametrize("raw_id", ACCEPTED_FORMS)
class TestAcceptedFormats:
    """Every spelling of one ID collapses to the same canonical forms."""

    def test_parses_to_the_number(self, raw_id):
        assert parse_chebi_number(raw_id) == "17079"

    def test_recognised_as_a_chebi_id(self, raw_id):
        assert looks_like_chebi_id(raw_id)

    def test_normalizes_to_one_iri(self, raw_id):
        assert normalize_id(raw_id) == f"{CHEBI_IRI_PREFIX}CHEBI_17079"

    def test_not_mistaken_for_smiles(self, raw_id):
        """The letters of 'chebi' are all SMILES atoms, so this is easy to get wrong."""
        assert not is_smiles(raw_id)


@pytest.mark.parametrize(
    "value",
    ["", "   ", "CHEBI:", "CHEBI ID:", "CCO", "CC(=O)Oc1ccccc1C(=O)O", "aspirin"],
)
def test_non_ids_are_not_parsed(value):
    assert parse_chebi_number(value) is None
    assert to_chebi_curie(value) is None
    assert to_chebi_iri(value) is None


def test_none_is_not_a_chebi_id():
    assert parse_chebi_number(None) is None


def test_non_chebi_iri_left_alone():
    """Only the ChEBI namespace is claimed; other IRIs pass through normalize_id."""
    iri = "http://purl.obolibrary.org/obo/GO_0008150"
    assert parse_chebi_number(iri) is None
    assert normalize_id(iri) == iri


def test_curie_form():
    assert to_chebi_curie("chebi id: 17079") == "CHEBI:17079"


def test_normalize_id_idempotent():
    once = normalize_id("chebi id: 17079")
    assert normalize_id(once) == once


class TestCollapseChebiPrefixes:
    """Prefixes containing whitespace have to survive the split into tokens."""

    def test_spaced_prefix_becomes_one_token(self):
        assert collapse_chebi_prefixes("CHEBI ID: 17079") == "CHEBI:17079"

    def test_keeps_the_rest_of_the_line(self):
        assert (
            collapse_chebi_prefixes("CHEBI ID: 17079\t0.7665") == "CHEBI:17079\t0.7665"
        )

    def test_several_ids_per_line(self):
        assert (
            collapse_chebi_prefixes("chebi 17079, chebi 17080")
            == "CHEBI:17079, CHEBI:17080"
        )

    def test_iris_are_left_intact(self):
        """Rewriting the CHEBI_ inside an IRI would break it."""
        iri = "http://purl.obolibrary.org/obo/CHEBI_17079"
        assert collapse_chebi_prefixes(iri) == iri

    def test_bare_numbers_untouched(self):
        assert collapse_chebi_prefixes("17079 17080") == "17079 17080"

    def test_smiles_untouched(self):
        smiles = "CC(=O)Oc1ccccc1C(=O)O 1.5"
        assert collapse_chebi_prefixes(smiles) == smiles
