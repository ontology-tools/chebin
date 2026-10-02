"""Tests for where the leaf InChIKeys come from.

The SMILES below are as ChEBI release 255 publishes them for uracil and
L-tryptophan. RDKit can't compute an InChIKey from either, so a SMILES or InChI
submitted for these compounds could only be matched through ChEBI's asserted key.
"""

import pandas as pd

from chebin.preparing_data.wikidata.get_inchikeys import (
    build_inchikey_map_from_owl,
    convert_smiles_file,
)

OBO = "http://purl.obolibrary.org/obo/"
URACIL_KEY = "ISAKRJDGNUQOIC-UHFFFAOYSA-N"
TRYPTOPHAN_KEY = "QIVBCDIJIAJPQS-VIFPVBQESA-N"
ETHANOL_KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

OWL = f"""<?xml version="1.0"?>
<rdf:RDF xmlns:owl="http://www.w3.org/2002/07/owl#"
         xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
         xmlns:chemrof="https://w3id.org/chemrof/">
    <owl:Class rdf:about="{OBO}CHEBI_17568">
        <rdfs:subClassOf rdf:resource="{OBO}CHEBI_26432"/>
        <chemrof:inchi_key_string>{URACIL_KEY}</chemrof:inchi_key_string>
        <chemrof:smiles_string>O=c1ccnc(=O)n1</chemrof:smiles_string>
    </owl:Class>
    <owl:Class rdf:about="{OBO}CHEBI_16828">
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:someValuesFrom rdf:resource="{OBO}CHEBI_1"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <chemrof:inchi_key_string>{TRYPTOPHAN_KEY}</chemrof:inchi_key_string>
        <chemrof:smiles_string>N[C@@H](Cc1cnc2ccccc12)C(=O)O</chemrof:smiles_string>
    </owl:Class>
    <owl:Class rdf:about="{OBO}CHEBI_16236">
        <chemrof:smiles_string>CCO</chemrof:smiles_string>
    </owl:Class>
</rdf:RDF>
"""


def _write_inputs(tmp_path):
    owl_file = tmp_path / "chebi.owl"
    owl_file.write_text(OWL, encoding="utf-8")
    leaves_csv = tmp_path / "leaves.csv"
    pd.DataFrame(
        {
            "IRI": [f"{OBO}CHEBI_17568", f"{OBO}CHEBI_16828", f"{OBO}CHEBI_16236"],
            "SMILES": [
                '"O=c1ccnc(=O)n1"',
                '"N[C@@H](Cc1cnc2ccccc12)C(=O)O"',
                '"CCO"',
            ],
            "Classification": ["structural"] * 3,
        },
    ).to_csv(leaves_csv, index=False)
    return owl_file, leaves_csv


def test_owl_map_reads_asserted_inchikeys(tmp_path):
    owl_file, _ = _write_inputs(tmp_path)
    assert build_inchikey_map_from_owl(owl_file) == {
        f"{OBO}CHEBI_17568": URACIL_KEY,
        f"{OBO}CHEBI_16828": TRYPTOPHAN_KEY,
    }


def test_asserted_inchikeys_are_used_and_computed_ones_fill_gaps(tmp_path):
    owl_file, leaves_csv = _write_inputs(tmp_path)
    out = tmp_path / "out.csv"
    convert_smiles_file(leaves_csv, out, owl_file=owl_file)
    keys = dict(zip(*pd.read_csv(out)[["IRI", "InChIKey"]].T.values))
    assert keys == {
        f"{OBO}CHEBI_17568": URACIL_KEY,
        f"{OBO}CHEBI_16828": TRYPTOPHAN_KEY,
        # No asserted key, so it is computed from the SMILES as before.
        f"{OBO}CHEBI_16236": ETHANOL_KEY,
    }


def test_without_owl_keys_are_only_computed(tmp_path):
    _, leaves_csv = _write_inputs(tmp_path)
    out = tmp_path / "out.csv"
    convert_smiles_file(leaves_csv, out)
    keys = pd.read_csv(out).set_index("IRI")["InChIKey"]
    assert keys.isna().sum() == 2
    assert keys[f"{OBO}CHEBI_16236"] == ETHANOL_KEY
