import pytest

from chebin.calculations import smiles_lookup


@pytest.fixture
def empty_local_maps(monkeypatch):
    """Stand in for the local SMILES/InChIKey -> ChEBI lookup tables with empty ones.

    Loading the real tables needs a generated data folder, which CI does not have.
    The structures in the tests that use this are either unreadable or absent from
    the real tables anyway, so they take the same path through
    convert_smiles_to_chebi with or without the data folder.
    """
    monkeypatch.setattr(smiles_lookup, "_local_maps", ({}, {}, {}, {}))
