import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import InchiToInchiKey, MolToInchi


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


def convert_smiles_file(input_file, output_file):

    # Open the input CSV file and read it
    df = pd.read_csv(input_file)
    print(f"Read {len(df)} rows from {input_file}")
    # Strip triple quotes from the SMILES column
    df["SMILES"] = df["SMILES"].str.strip('"')
    # Canonicalize (so the website's exact-string lookup is toolkit-consistent)
    # and convert to InChIKey in one pass per row.
    converted = df["SMILES"].apply(_canonicalize_and_convert)
    df["SMILES"] = converted.apply(lambda t: t[0])
    df["InChIKey"] = converted.apply(lambda t: t[1])
    starcount = converted.apply(lambda t: t[2]).sum()
    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_file, index=False)

    # Print summary statistics
    total_rows = len(df)
    generated_keys = df["InChIKey"].notnull().sum()
    failed_conversions = total_rows - generated_keys
    print(f"Processed {total_rows} rows.")
    print(f"Generated {generated_keys} InChIKeys.")
    print(f"Failed conversions: {failed_conversions}.")
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

    convert_smiles_file(input_file, output_file)
    count_nans(output_file)
