import argparse
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import cast

import pandas as pd
import requests

from chebin.calculations.smiles_lookup import CHEBIFIER_CLASSIFY_MODEL

# convert_smiles_to_chebi is a blocking HTTP call to a third-party API, so rows
# resolve concurrently on a thread pool rather than one at a time. Kept modest
# to stay a well-behaved client of a service we don't control.
_RESOLVE_WORKERS = 8


def normalize_chebi_id(raw_value):
    if pd.isna(raw_value):
        return None
    value = str(raw_value).strip()
    if not value:
        return None

    normalized = []
    for part in re.split(r"[|;]", value):
        part = part.strip()
        if part.startswith("http://purl.obolibrary.org/obo/CHEBI_"):
            normalized.append(part.rsplit("/", 1)[-1])
        elif part.startswith("CHEBI:"):
            normalized.append(part.replace(":", "_", 1))
        elif part.startswith("CHEBI_"):
            normalized.append(part)
        elif part.isdigit():
            normalized.append(f"CHEBI_{part}")
        else:
            normalized.append(part)

    return "|".join(normalized)


def convert_smiles_to_chebi(smiles_string):

    chebi_ids = []
    was_resolved_directly = False
    parents_found = False

    # Get details from ChEBI lookup to check for a direct match to a ChEBI ID.
    response = requests.post(
        "https://chebifier.hastingslab.org/api/details",
        json={
            "type": "type",
            "smiles": smiles_string,
            "selectedModels": {
                "ChEBI Lookup": True,
            },
        },
    )

    lookup_model = response.json().get("models", {}).get("ChEBI Lookup", {})

    # Prefer the API's structured chebi_ids list, falling back to pulling the
    # IDs out of the human-readable highlights text.
    lookup_ids = lookup_model.get("chebi_ids")
    if not lookup_ids:
        lookup_infotext = lookup_model.get("highlights", [])
        if lookup_infotext:
            lookup_ids = re.findall(r"CHEBI:(\d+)", lookup_infotext[0][1])

    # If the lookup found a ChEBI ID, use that directly. A SMILES matching
    # several entries contributes only one of them, tie-broken on the lowest
    # ChEBI ID to match the website's behaviour.
    if lookup_ids:
        chosen_id = min(lookup_ids, key=int)
        chebi_ids.append(f"CHEBI:{chosen_id}")
        was_resolved_directly = True
        # print(f"Found ChEBI ID from lookup: CHEBI:{chosen_id} for SMILES {smiles_string}")
    else:
        # print(f"No direct ChEBI ID found from lookup for SMILES {smiles_string}, attempting classification...")
        response = requests.post(
            "https://chebifier.hastingslab.org/api/classify",
            json={
                "smiles": smiles_string,
                "ontology": False,
                "selectedModels": {
                    CHEBIFIER_CLASSIFY_MODEL: True,
                },
            },
        )

        direct_parents = response.json().get("direct_parents")
        if direct_parents:
            for parent_list in direct_parents:
                if parent_list is not None:
                    parent_ids = [f"CHEBI:{parent[0]}" for parent in parent_list]
                    chebi_ids.extend(parent_ids)
                    parents_found = True
                    # print(f"Found direct parent ChEBI IDs from classification for SMILES {smiles_string}: {parent_ids}")
                else:
                    # print(f"No parents found in one of the classification results for SMILES {smiles_string}")
                    # print(f"Classification response content: {response.content}")
                    pass

            if chebi_ids:
                # print(f"Found {len(chebi_ids)} ChEBI IDs from classification for SMILES {smiles_string}")
                pass

    if not chebi_ids:
        # print(f"No ChEBI IDs found for SMILES {smiles_string} after lookup and classification.")
        pass

    return chebi_ids, was_resolved_directly, parents_found


def pick_smiles_candidates(raw_value):
    if pd.isna(raw_value):
        return []

    value = str(raw_value).strip()
    if not value:
        return []

    # Some rows contain multiple alternative SMILES joined by '|'.
    # Return all non-empty entries in order so we can try them one by one.
    candidates = []
    for part in value.split("|"):
        part = part.strip()
        if part:
            candidates.append(part)

    return candidates


def _choose_smiles_columns(df, requested_smiles_columns=None):
    """Pick SMILES columns that exist in the file.

    If `requested_smiles_columns` is provided, only those are considered and
    must exist. Otherwise, auto-detect common schemas.
    """

    if requested_smiles_columns:
        missing = [col for col in requested_smiles_columns if col not in df.columns]
        if missing:
            raise KeyError(f"Requested SMILES columns not found: {missing}")
        return list(requested_smiles_columns)

    auto_candidates = ["canonicalSmiles", "isomericSmiles", "smiles"]
    chosen = [col for col in auto_candidates if col in df.columns]
    if not chosen:
        raise KeyError(
            "Could not find a SMILES column. Expected at least one of "
            "['canonicalSmiles', 'isomericSmiles', 'smiles'].",
        )
    return chosen


def _resolve_smiles_candidates(smiles_candidates):
    """Resolve one row's SMILES candidates to ChEBI IDs.

    Tries each candidate in order, preferring a direct match; falls back to
    the first candidate's parent-based match if no candidate resolves
    directly. Runs on the thread pool in find_missing_chebis, so it must not
    touch the DataFrame -- only the network and its own locals.
    """
    chebi_ids = []
    was_resolved_directly = False
    parents_found = False
    first_parent_ids = None

    for smiles_to_query in smiles_candidates:
        candidate_ids, candidate_direct, candidate_parents = convert_smiles_to_chebi(
            smiles_to_query,
        )

        if candidate_direct and candidate_ids:
            return candidate_ids, True, False

        if candidate_parents and candidate_ids and first_parent_ids is None:
            first_parent_ids = candidate_ids

    if first_parent_ids:
        chebi_ids = first_parent_ids
        parents_found = True

    return chebi_ids, was_resolved_directly, parents_found


def find_missing_chebis(
    compounds_file,
    output_file_path=None,
    smiles_columns=None,
    chebi_column="chebi_id",
):
    df = pd.read_csv(compounds_file, sep="\t")

    if chebi_column not in df.columns:
        raise KeyError(f"Missing required column '{chebi_column}' in {compounds_file}")

    selected_smiles_columns = _choose_smiles_columns(df, smiles_columns)
    print(f"Using SMILES columns: {selected_smiles_columns}")

    if "chebi_source" not in df.columns:
        df["chebi_source"] = None

    # Standardize existing IDs to CHEBI_##### format.
    df[chebi_column] = df[chebi_column].apply(normalize_chebi_id)

    # Treat both NaN and empty strings as missing ChEBI IDs.
    chebi_existing_mask = df[chebi_column].notna() & (
        df[chebi_column].astype(str).str.strip() != ""
    )
    df.loc[chebi_existing_mask, "chebi_source"] = "already_existed"

    missing_indices = df.index[~chebi_existing_mask]
    print(f"Rows with missing ChEBI IDs: {len(missing_indices)}")

    # Build optional InChIKey -> ChEBI IRI map from removed leaf classes file.
    inchikey_map = {}
    try:
        removed_inchikeys = pd.read_csv("data/removed_leaf_classes_with_inchikeys.csv")
        # possible column names for InChIKey: InChIKey, InChIkey, inchikey
        inchikey_col = None
        for c in ["InChIKey", "InChIkey", "inchikey", "InChIKEY"]:
            if c in removed_inchikeys.columns:
                inchikey_col = c
                break
        if inchikey_col is not None and "IRI" in removed_inchikeys.columns:
            for _, r in removed_inchikeys.dropna(subset=[inchikey_col]).iterrows():
                key = str(r[inchikey_col]).strip()
                if key:
                    # map InChIKey -> normalized CHEBI id (from IRI)
                    iri = r.get("IRI")
                    if pd.notna(iri) and iri:
                        inchikey_map[key] = normalize_chebi_id(iri)
        if inchikey_map:
            print(
                f"Loaded {len(inchikey_map)} InChIKey -> ChEBI mappings from removed_leaf_classes_with_inchikeys.csv",
            )
    except FileNotFoundError:
        pass

    # Counter for how many rows were resolved via InChIKey mapping.
    inchikey_match_count = 0

    # First pass (cheap, local only): resolve what we can from the InChIKey
    # map, and work out which rows actually need an API call.
    pending_candidates = {}
    for idx in missing_indices:
        smiles_candidates = []
        for smiles_column in selected_smiles_columns:
            candidates = pick_smiles_candidates(df.at[idx, smiles_column])
            for candidate in candidates:
                if candidate not in smiles_candidates:
                    smiles_candidates.append(candidate)

        # If an InChIKey column exists in the input and maps to a ChEBI, use it first.
        inchikey_used = False
        for colname in ("InChIKey", "InChIkey", "inchikey", "inchiKey"):
            if colname in df.columns:
                ik = str(df.at[idx, colname]).strip()
                if ik and ik in inchikey_map:
                    df.at[idx, chebi_column] = inchikey_map[ik]
                    df.at[idx, "chebi_source"] = "found_via_inchikey"
                    inchikey_match_count += 1
                    inchikey_used = True
                break

        if inchikey_used:
            continue

        if not smiles_candidates:
            df.at[idx, "chebi_source"] = "unresolved"
            continue

        pending_candidates[idx] = smiles_candidates

    # Second pass: resolve the remaining rows concurrently. Each row is an
    # independent 1-2 request round trip to the Chebifier API, so this is
    # waiting on network latency rather than doing CPU work -- a thread pool
    # lets those round trips overlap instead of running one at a time.
    print(f"Resolving {len(pending_candidates)} rows via Chebifier API...")
    processed = 0
    with ThreadPoolExecutor(max_workers=_RESOLVE_WORKERS) as executor:
        future_to_idx = {
            executor.submit(_resolve_smiles_candidates, candidates): idx
            for idx, candidates in pending_candidates.items()
        }
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            chebi_ids, was_resolved_directly, parents_found = future.result()

            if chebi_ids:
                df.at[idx, chebi_column] = normalize_chebi_id("|".join(chebi_ids))
                if was_resolved_directly:
                    df.at[idx, "chebi_source"] = "found_directly"
                elif parents_found:
                    df.at[idx, "chebi_source"] = "parents_compounds"
                else:
                    df.at[idx, "chebi_source"] = "unresolved"
            else:
                df.at[idx, "chebi_source"] = "unresolved"

            processed += 1
            if processed % 50 == 0:
                print(f"Processed {processed}/{len(pending_candidates)} missing rows")

    save_path = output_file_path if output_file_path else compounds_file
    df.to_csv(save_path, sep="\t", index=False)
    print(f"Saved updated compounds file to {save_path}")

    direct_count = int((df["chebi_source"] == "found_directly").sum())
    parents_count = int((df["chebi_source"] == "parents_compounds").sum())
    still_missing_mask = df[chebi_column].isna() | (
        df[chebi_column].astype(str).str.strip() == ""
    )
    unresolved_count = int(still_missing_mask.sum())

    print(f"Resolved via InChIKey mapping: {inchikey_match_count}")

    print(f"Found directly: {direct_count}")
    print(f"Found using parents: {parents_count}")
    print(f"Still unresolved (no chebi_id): {unresolved_count}")

    return df


def run_find_missing_chebis(
    source: str = "wikidata_hs",
    compounds_file: str | None = None,
    output_file: str | None = None,
    smiles_columns: list | None = None,
    chebi_column: str | None = None,
):
    """Programmatic wrapper around :func:`find_missing_chebis`.

    Parameters mirror the CLI presets. This avoids argparse when calling
    from other scripts (e.g., `create_files.py`).
    """
    SOURCE_PRESETS = {
        "wikidata_hs": {
            "input": "data/intermediate_files/compounds_with_chebi_ids_homo_sapiens.tsv",
            "output": "data/intermediate_files/compounds_with_chebi_ids_homo_sapiens_updatedchebis.tsv",
            "smiles_columns": ["canonicalSmiles", "isomericSmiles"],
            "chebi_column": "chebi_id",
        },
        "hmdb": {
            "input": "data/hmdb_metabolites_extract_quantified_detected.tsv",
            "output": "data/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv",
            "smiles_columns": ["smiles"],
            "chebi_column": "chebi_id",
        },
        "wikidata_at": {
            "input": "data/intermediate_files/compounds_with_chebi_ids_arabidopsis_thaliana.tsv",
            "output": "data/intermediate_files/compounds_with_chebi_ids_arabidopsis_thaliana_updatedchebis.tsv",
            "smiles_columns": ["canonicalSmiles", "isomericSmiles"],
            "chebi_column": "chebi_id",
        },
        "lotus_hs": {
            "input": "data/intermediate_files/lotus_homo_sapiens_with_chebi_ids.tsv",
            "output": "data/intermediate_files/lotus_homo_sapiens_with_chebi_ids_updatedchebis.tsv",
            "smiles_columns": ["compound_smiles_conn", "compound_smiles_iso"],
            "chebi_column": "chebi_id",
        },
        "lotus_at": {
            "input": "data/intermediate_files/lotus_arabidopsis_thaliana_with_chebi_ids.tsv",
            "output": "data/intermediate_files/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv",
            "smiles_columns": ["compound_smiles_conn", "compound_smiles_iso"],
            "chebi_column": "chebi_id",
        },
    }

    if source not in SOURCE_PRESETS:
        raise ValueError(f"Unknown source preset: {source!r}")

    preset = SOURCE_PRESETS[source]
    compounds_file = (
        compounds_file if compounds_file is not None else cast(str, preset["input"])
    )
    output_file = (
        output_file if output_file is not None else cast(str, preset["output"])
    )
    smiles_columns = (
        smiles_columns
        if smiles_columns is not None
        else cast(list[str], preset["smiles_columns"])
    )
    chebi_column = (
        chebi_column if chebi_column is not None else cast(str, preset["chebi_column"])
    )

    print(
        f"Running find_missing_chebis programmatically: source={source}, input={compounds_file}, output={output_file}",
    )
    return find_missing_chebis(
        compounds_file,
        output_file,
        smiles_columns=smiles_columns,
        chebi_column=chebi_column,
    )


def main_find_missing_chebis(source):
    start_time = time.time()

    SOURCE_PRESETS = {
        "wikidata_hs": {
            "input": "data/intermediate_files/compounds_with_chebi_ids_homo_sapiens.tsv",
            "output": "data/intermediate_files/compounds_with_chebi_ids_homo_sapiens_updatedchebis.tsv",
            "smiles_columns": ["canonicalSmiles", "isomericSmiles"],
            "chebi_column": "chebi_id",
        },
        "hmdb": {
            "input": "data/hmdb_metabolites_extract_quantified_detected.tsv",
            "output": "data/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv",
            "smiles_columns": ["smiles"],
            "chebi_column": "chebi_id",
        },
        "wikidata_at": {
            "input": "data/intermediate_files/compounds_with_chebi_ids_arabidopsis_thaliana.tsv",
            "output": "data/intermediate_files/compounds_with_chebi_ids_arabidopsis_thaliana_updatedchebis.tsv",
            "smiles_columns": ["canonicalSmiles", "isomericSmiles"],
            "chebi_column": "chebi_id",
        },
        "lotus_hs": {
            "input": "data/intermediate_files/lotus_homo_sapiens_with_chebi_ids.tsv",
            "output": "data/intermediate_files/lotus_homo_sapiens_with_chebi_ids_updatedchebis.tsv",
            "smiles_columns": ["compound_smiles_conn", "compound_smiles_iso"],
            "chebi_column": "chebi_id",
        },
        "lotus_at": {
            "input": "data/intermediate_files/lotus_arabidopsis_thaliana_with_chebi_ids.tsv",
            "output": "data/intermediate_files/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv",
            "smiles_columns": ["compound_smiles_conn", "compound_smiles_iso"],
            "chebi_column": "chebi_id",
        },
    }

    parser = argparse.ArgumentParser(
        description="Fill missing ChEBI IDs using SMILES for Wikidata/HMDB-style TSV files.",
    )
    parser.add_argument(
        "--source",
        choices=("wikidata_hs", "hmdb", "wikidata_at", "lotus_hs", "lotus_at"),
        default=source,
        help="Use built-in defaults for selected source (default: wikidata_hs).",
    )
    parser.add_argument(
        "compounds_file",
        nargs="?",
        default=None,
        help="Input TSV file with at least a chebi_id column and one SMILES column.",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        default=None,
        help="Output TSV file path.",
    )
    parser.add_argument(
        "--smiles-columns",
        nargs="+",
        default=None,
        help="Optional explicit SMILES columns (e.g., --smiles-columns smiles or canonicalSmiles isomericSmiles).",
    )
    parser.add_argument(
        "--chebi-column",
        default=None,
        help="Name of the ChEBI ID column (default: chebi_id).",
    )

    args = parser.parse_args()

    preset = SOURCE_PRESETS[args.source]
    compounds_file = args.compounds_file if args.compounds_file else preset["input"]
    output_file = args.output_file if args.output_file else preset["output"]
    smiles_columns = (
        args.smiles_columns if args.smiles_columns else preset["smiles_columns"]
    )
    chebi_column = args.chebi_column if args.chebi_column else preset["chebi_column"]

    print(f"Source preset: {args.source}")
    print(f"Input file: {compounds_file}")
    print(f"Output file: {output_file}")

    find_missing_chebis(
        compounds_file,
        output_file,
        smiles_columns=smiles_columns,
        chebi_column=chebi_column,
    )

    elapsed_seconds = time.time() - start_time
    elapsed_minutes = elapsed_seconds / 60
    print(
        f"Total runtime: {elapsed_seconds:.2f} seconds ({elapsed_minutes:.2f} minutes)",
    )


if __name__ == "__main__":
    source = "wikidata_hs"  # or "hmdb" or "wikidata_at"
    main_find_missing_chebis(source)
