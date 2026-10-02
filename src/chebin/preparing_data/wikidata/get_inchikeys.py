import xml.etree.ElementTree as ET

import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi

_OWL_CLASS_TAG = "{http://www.w3.org/2002/07/owl#}Class"
_RDF_ABOUT = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
_INCHIKEY_TAG = "{https://w3id.org/chemrof/}inchi_key_string"


def build_inchikey_map_from_owl(owl_file):
    """Map each ChEBI class IRI to the InChIKey ChEBI asserts for it.

    ChEBI's asserted InChIKey is preferred over one recomputed from the class's
    SMILES with RDKit, since the two don't always agree -- and for some SMILES
    RDKit can't compute one at all.
    """
    inchikey_map = {}
    for _, elem in ET.iterparse(owl_file, events=("end",)):
        if elem.tag == _OWL_CLASS_TAG:
            iri = elem.get(_RDF_ABOUT)
            key_elem = elem.find(_INCHIKEY_TAG)
            if iri and key_elem is not None and key_elem.text:
                inchikey_map[iri] = key_elem.text.strip()
            elem.clear()
    return inchikey_map


def smiles_to_inchikey(smiles):
    """InChIKey for a SMILES string, or None if it can't be parsed or has a star atom.

    Returns (inchikey_or_None, had_star) so callers can tally star counts
    without a second RDKit parse.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, False
    if "*" in smiles:
        return None, True
    inchi = MolToInchi(mol)
    return InchiToInchiKey(inchi), False


def canonical_smiles(smiles):
    """RDKit-canonical form of a SMILES string, or None if it can't be parsed.

    ChEBI's asserted SMILES aren't guaranteed to be in any particular canonical
    form, so baking the RDKit-canonical form into this file lets the website's
    exact-string lookup match incoming SMILES regardless of how they were written.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def _canonicalize_and_convert(smiles):
    """Canonical SMILES, InChIKey, and star-atom flag from a single RDKit parse.

    Canonicalizing and converting separately (canonical_smiles then
    smiles_to_inchikey) parses the same molecule with RDKit twice; doing both
    off one parsed Mol roughly triples throughput on large files.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, None, False
    canonical = Chem.MolToSmiles(mol)
    if "*" in smiles:
        return canonical, None, True
    return canonical, InchiToInchiKey(MolToInchi(mol)), False


def convert_smiles_file(input_file, output_file, owl_file=None):
    """Add canonical SMILES and an InChIKey to the removed-leaf-classes CSV.

    owl_file: the ChEBI OWL the leaves came from. When given, each leaf gets the
    InChIKey ChEBI asserts for it (see build_inchikey_map_from_owl), and one is
    only computed from the SMILES for leaves without an asserted InChIKey.
    """

    # Open the input CSV file and read it
    df = pd.read_csv(input_file)
    print(f"Read {len(df)} rows from {input_file}")
    # Strip triple quotes from the SMILES column
    df["SMILES"] = df["SMILES"].str.strip('"')
    # Canonicalize (so the website's exact-string lookup is toolkit-consistent)
    # and convert to InChIKey in one pass per row.
    converted = df["SMILES"].apply(_canonicalize_and_convert)
    df["SMILES"] = converted.apply(lambda t: t[0])
    computed_keys = converted.apply(lambda t: t[1])
    starcount = converted.apply(lambda t: t[2]).sum()

    if owl_file is not None:
        owl_keys = df["IRI"].map(build_inchikey_map_from_owl(owl_file))
        df["InChIKey"] = owl_keys.fillna(computed_keys)
        both = owl_keys.notna() & computed_keys.notna()
        print(f"InChIKeys asserted by ChEBI: {int(owl_keys.notna().sum())}.")
        print(
            f"  of which none could be computed from the SMILES: "
            f"{int((owl_keys.notna() & computed_keys.isna()).sum())}, "
            f"and the computed one differed: "
            f"{int((owl_keys[both] != computed_keys[both]).sum())}.",
        )
    else:
        df["InChIKey"] = computed_keys
    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_file, index=False)

    # Print summary statistics
    total_rows = len(df)
    generated_keys = df["InChIKey"].notnull().sum()
    failed_conversions = total_rows - generated_keys
    print(f"Processed {total_rows} rows.")
    print(f"Rows with an InChIKey: {generated_keys}.")
    print(f"Rows without one: {failed_conversions}.")
    print(f"SMILES with stars (not converted): {starcount}.")


def count_nans(csv_file):
    df = pd.read_csv(csv_file)
    # check how many of the rows that do not have an inchikey that have Classification 'strictural'
    column_name = "InChIKey"
    na_count = df[column_name].isna().sum()
    structural_na_count = (
        df[df["Classification"] == "structural"][column_name].isna().sum()
    )

    print(f"Total NaN in '{column_name}': {na_count}")
    print(
        f"NaN in '{column_name}' where Classification is 'structural': {structural_na_count}",
    )


if __name__ == "__main__":
    input_file = "data/removed_leaf_classes_with_smiles.csv"
    output_file = "data/removed_leaf_classes_with_inchikeys.csv"

    convert_smiles_file(input_file, output_file, owl_file="data/source_files/chebi.owl")
    count_nans(output_file)
