"""find_missing_chebis must read the InChIKey table it is given.

create_all_files_with_backup builds into a temporary folder, so a hardcoded
``data/...`` path would silently match against the previous build's table.
"""

import pandas as pd
import pytest

from chebin.preparing_data.wikidata import find_missing_chebis as fmc

URACIL_KEY = "ISAKRJDGNUQOIC-UHFFFAOYSA-N"


@pytest.fixture
def no_network(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("made a network call")

    monkeypatch.setattr(fmc.requests, "post", forbidden)


def test_reads_the_given_inchikey_table(tmp_path, monkeypatch, no_network):
    # Run from an empty folder, so there is no data/ table to fall back on.
    monkeypatch.chdir(tmp_path)
    inchikeys_csv = tmp_path / "new_build" / "removed_leaf_classes_with_inchikeys.csv"
    inchikeys_csv.parent.mkdir()
    pd.DataFrame(
        {
            "IRI": ["http://purl.obolibrary.org/obo/CHEBI_17568"],
            "SMILES": ["O=c1ccnc(=O)n1"],
            "InChIKey": [URACIL_KEY],
        },
    ).to_csv(inchikeys_csv, index=False)
    compounds = tmp_path / "compounds.tsv"
    pd.DataFrame(
        {
            "chebi_id": [None],
            "smiles": ["O=c1cc[nH]c(=O)[nH]1"],
            "InChIKey": [URACIL_KEY],
        },
    ).to_csv(compounds, sep="\t", index=False)

    df = fmc.find_missing_chebis(
        compounds,
        tmp_path / "out.tsv",
        smiles_columns=["smiles"],
        inchikeys_csv=str(inchikeys_csv),
    )

    assert df.loc[0, "chebi_id"] == "CHEBI_17568"
    assert df.loc[0, "chebi_source"] == "found_via_inchikey"
