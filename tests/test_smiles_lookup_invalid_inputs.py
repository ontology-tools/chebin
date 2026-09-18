"""Tests for how smiles_lookup handles input it cannot resolve.

Two failures used to look the same from outside and both ended badly: an input
RDKit can't read was still sent to Chebifier, which answers some of them with a
502 HTML page, and the unguarded ``.json()`` on that page raised straight through
to the user as a 500. Neither the wasted call nor the crash should happen, and the
two causes have to stay distinguishable -- "your input is a typo" and "the lookup
service is down" need different answers from the user.

Every test here stubs ``requests.post``; nothing in this file touches the network.
"""

import pytest
import requests

from chebin.calculations import smiles_lookup

# Unreadable by RDKit. The exact string from the bug report.
INVALID_SMILES = "CCO=OOOO"
# Valid, and absent from the local lookup tables, so resolving it needs the API.
VALID_REMOTE_SMILES = "CCOC(=O)CCCCCCN1CCCC1"


@pytest.fixture
def no_network(monkeypatch):
    """Fail loudly if anything reaches the network."""

    def forbidden(*args, **kwargs):
        raise AssertionError("made a network call")

    monkeypatch.setattr(smiles_lookup.requests, "post", forbidden)


def _stub_response(monkeypatch, *, status_error=None, json_error=None):
    class Response:
        content = b"<html><head><title>502 Bad Gateway</title></head></html>"

        def raise_for_status(self):
            if status_error is not None:
                raise status_error

        def json(self):
            if json_error is not None:
                raise json_error
            return {}

    monkeypatch.setattr(smiles_lookup.requests, "post", lambda *a, **k: Response())


@pytest.mark.parametrize("use_parents", [False, True])
def test_unreadable_input_never_reaches_the_api(no_network, use_parents):
    """RDKit already knows it isn't a molecule, so there is nothing to ask about."""
    chebi_ids, was_resolved, ambiguous, failure_reason = (
        smiles_lookup.convert_smiles_to_chebi(
            INVALID_SMILES,
            use_parents=use_parents,
        )
    )
    assert chebi_ids == []
    assert was_resolved is False
    assert ambiguous is None
    assert failure_reason == "invalid_structure"


def test_unreadable_inchi_is_also_invalid(no_network):
    _, was_resolved, _, failure_reason = smiles_lookup.convert_smiles_to_chebi(
        "InChI=1S/not-a-real-inchi",
    )
    assert was_resolved is False
    assert failure_reason == "invalid_structure"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"status_error": requests.HTTPError("502 Server Error")},
        {"json_error": requests.exceptions.JSONDecodeError("Expecting value", "", 0)},
        {"status_error": requests.Timeout("read timed out")},
    ],
    ids=["http-502", "html-body-instead-of-json", "timeout"],
)
def test_api_failure_degrades_instead_of_raising(monkeypatch, kwargs):
    """A broken lookup must not be reported as invalid input -- the structure is fine."""
    _stub_response(monkeypatch, **kwargs)
    chebi_ids, was_resolved, _, failure_reason = smiles_lookup.convert_smiles_to_chebi(
        VALID_REMOTE_SMILES,
        use_parents=True,
    )
    assert chebi_ids == []
    assert was_resolved is False
    assert failure_reason == "lookup_unavailable"


def test_wrappers_report_invalid_as_a_subset_of_unresolved(no_network):
    """Callers reading only unresolved_smiles must still see every failed input."""
    studyset, unresolved, ambiguous, invalid = smiles_lookup.smiles_list_to_studyset(
        [INVALID_SMILES, "CHEBI:17079"],
    )
    assert studyset == ["CHEBI:17079"]
    assert unresolved == [INVALID_SMILES]
    assert invalid == [INVALID_SMILES]
    assert ambiguous == []


def test_weights_wrapper_reports_invalid_structures(no_network):
    weights, unresolved, _, invalid = smiles_lookup.smiles_weights_to_chebi_weights(
        {INVALID_SMILES: 0.5, "CHEBI:17079": 2.0},
    )
    assert weights == {"CHEBI:17079": 2.0}
    assert unresolved == [INVALID_SMILES]
    assert invalid == [INVALID_SMILES]
