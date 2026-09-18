"""Tests for the website's study-set parser (website/app.py parse_studyset).

This is where a pasted study set turns into the ID list the analysis runs on, so
it is where a ChEBI ID written in an unexpected format used to be lost: the whole
run came back empty with nothing to say why. The parser is also the only place
that has to deal with IDs whose prefix contains a space ("CHEBI ID: 17079"),
since it splits lines on whitespace to pick up the optional weight column.
"""

import os
import sys
from pathlib import Path

import pytest

WEBSITE_DIR = Path(__file__).resolve().parents[1] / "website"


@pytest.fixture(scope="module")
def website_app():
    """Import website/app.py, leaving the working directory as we found it.

    app.py chdirs to the project root and points chebin at its data folder on
    import (the website relies on both), neither of which should leak into the
    rest of the suite.
    """
    cwd = os.getcwd()
    sys.path.insert(0, str(WEBSITE_DIR))
    try:
        import app
    finally:
        sys.path.remove(str(WEBSITE_DIR))
        os.chdir(cwd)
    return app


@pytest.fixture
def parse(website_app):
    """parse_studyset() inside a request context, since it reads the session."""

    def _parse(text):
        with website_app.app.test_request_context():
            from flask import session

            # Keep SMILES resolution off the network: nothing here is a SMILES.
            session["smiles_option"] = "exclude"
            return website_app.parse_studyset(text)

    return _parse


@pytest.mark.parametrize(
    "text",
    [
        "CHEBI:17079",
        "chebi:17079",
        "ChEBI:17079",
        "CHEBI_17079",
        "CHEBI 17079",
        "CHEBI ID: 17079",
        "CHEBI ID:17079",
        "chebi id 17079",
        "17079",
        "  CHEBI:17079  ",
        "http://purl.obolibrary.org/obo/CHEBI_17079",
    ],
)
def test_every_id_format_reaches_the_study_set(parse, text):
    studyset, weights, unresolved, ambiguous = parse(text)
    assert studyset == ["CHEBI_17079"]
    assert weights == {}
    assert unresolved == []
    assert ambiguous == []


def test_mixed_formats_in_one_submission(parse):
    studyset, _, _, _ = parse("CHEBI:17079\nchebi:46816\nCHEBI ID: 31463\n28426")
    assert studyset == ["CHEBI_17079", "CHEBI_46816", "CHEBI_31463", "CHEBI_28426"]


def test_weights_survive_a_spaced_prefix(parse):
    """The weight column still lines up once "CHEBI ID: 17079" is one token."""
    studyset, weights, _, _ = parse("CHEBI ID: 17079\t0.7665\nchebi:46816 0.7465")
    assert studyset == ["CHEBI_17079", "CHEBI_46816"]
    assert weights == {"CHEBI_17079": 0.7665, "CHEBI_46816": 0.7465}


def test_bare_ids_with_weights(parse):
    studyset, weights, _, _ = parse("17079 0.7665")
    assert studyset == ["CHEBI_17079"]
    assert weights == {"CHEBI_17079": 0.7665}


def test_two_bare_ids_are_not_an_id_and_a_weight(parse):
    """Ambiguous on its own, and read as IDs: a weight has to be a decimal here."""
    studyset, weights, _, _ = parse("17079 17080")
    assert studyset == ["CHEBI_17079", "CHEBI_17080"]
    assert weights == {}


def test_comma_separated_bare_ids(parse):
    studyset, weights, _, _ = parse("17079, 17080, 17081")
    assert studyset == ["CHEBI_17079", "CHEBI_17080", "CHEBI_17081"]
    assert weights == {}


def test_prefixed_id_keeps_a_whole_number_weight(parse):
    """With a prefix there is no ambiguity, so an integer weight still counts."""
    studyset, weights, _, _ = parse("CHEBI:17079 2")
    assert studyset == ["CHEBI_17079"]
    assert weights == {"CHEBI_17079": 2.0}


def test_blank_input(parse):
    assert parse("   \n\n  ") == ([], {}, [], [])


GLUCOSE_INCHI = (
    "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1"
)


def test_inchi_is_one_entry(website_app):
    """An InChI contains commas, so a plain comma split scattered it into six."""
    assert website_app._split_entries(GLUCOSE_INCHI) == [GLUCOSE_INCHI]


def test_inchi_keeps_its_weight_column(website_app):
    parts = website_app._split_entries(f"{GLUCOSE_INCHI}\t0.5")
    assert parts == [GLUCOSE_INCHI, "0.5"]
    assert website_app._line_weight(parts) == 0.5


def test_commas_still_separate_non_inchi_entries(website_app):
    assert website_app._split_entries("CCO, CHEBI:17079,17080") == [
        "CCO",
        "CHEBI:17079",
        "17080",
    ]
