"""SMILES -> ChEBI ID resolution.

Shared by the ``_from_smiles`` enrichment wrappers (:mod:`chebin.calculations.fishers_calculations`,
:mod:`chebin.calculations.weighted_calculations`,
:mod:`chebin.preparing_data.wikidata.narrow_background_fishers`) and the website, so there's
one implementation of the lookup cascade instead of one per caller.
"""

from __future__ import annotations

import csv
import re

import requests
from rdkit import Chem
from rdkit.Chem import inchi

from chebin.calculations.chebi_ids import looks_like_chebi_id, to_chebi_curie
from chebin.calculations.log_utils import preview
from chebin.config import require_data_path

#: Local table asserting SMILES/InChIKey -> ChEBI ID for removed leaf classes, relative
#: to the data folder (see chebin.config).
SMILES_INCHIKEY_LOOKUP_CSV = "removed_leaf_classes_with_inchikeys.csv"

CHEBIFIER_DETAILS_URL = "https://chebifier.hastingslab.org/api/details"
CHEBIFIER_CLASSIFY_URL = "https://chebifier.hastingslab.org/api/classify"


def _clean_lookup_value(value):
    if value is None:
        return ""
    return str(value).strip().strip('"')


def _normalize_chebi_id(raw_value):
    value = _clean_lookup_value(raw_value)
    if not value:
        return ""

    curie = to_chebi_curie(value)
    if curie is not None:
        return curie.replace(":", "_", 1)

    return value


def _chebi_sort_key(chebi_id):
    """Numeric sort key for 'CHEBI_12345' ids, so collision tie-breaks are by
    ChEBI ID number rather than CSV row order (which is incidental).
    """
    try:
        return int(chebi_id.rsplit("_", 1)[-1])
    except (ValueError, AttributeError):
        return float("inf")


def _resolve_lookup_collisions(candidates, label):
    """Pick a deterministic winner (lowest ChEBI ID) for each key asserted by more
    than one ChEBI term, and return both the winner map and a {key: [all_ids]} map
    of only the keys that actually collided, so callers can surface the ones that
    were dropped. Logs a single summary line rather than one line per collision,
    since this table has thousands of them; per-match ambiguity is surfaced to end
    users via the ambiguous_matches returned by smiles_list_to_studyset /
    smiles_weights_to_chebi_weights instead.
    """
    resolved = {}
    collisions = {}
    for key, chebi_ids in candidates.items():
        ordered = sorted(chebi_ids, key=_chebi_sort_key)
        resolved[key] = ordered[0]
        if len(ordered) > 1:
            collisions[key] = ordered
    if collisions:
        print(
            f"Warning: {len(collisions)} distinct {label} values in {SMILES_INCHIKEY_LOOKUP_CSV} "
            f"are each asserted by more than one ChEBI term; using the lowest ChEBI ID "
            f"as a deterministic tie-break for each.",
        )
    return resolved, collisions


def _load_local_smiles_and_inchikey_maps():
    smiles_candidates = {}
    inchikey_candidates = {}

    lookup_file = require_data_path(SMILES_INCHIKEY_LOOKUP_CSV)
    with open(lookup_file, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []

        smiles_column = None
        inchikey_column = None
        iri_column = None

        for candidate in ("SMILES", "smiles"):
            if candidate in fieldnames:
                smiles_column = candidate
                break

        for candidate in ("InChIKey", "InChIkey", "inchikey", "InChIKEY"):
            if candidate in fieldnames:
                inchikey_column = candidate
                break

        for candidate in ("IRI", "iri"):
            if candidate in fieldnames:
                iri_column = candidate
                break

        if smiles_column is None or iri_column is None:
            raise KeyError(
                f"{lookup_file} must contain SMILES and IRI columns to build local lookup maps",
            )

        for row in reader:
            chebi_id = _normalize_chebi_id(row.get(iri_column))
            smiles = _clean_lookup_value(row.get(smiles_column))
            inchikey_value = (
                _clean_lookup_value(row.get(inchikey_column)) if inchikey_column else ""
            )

            if smiles:
                ids = smiles_candidates.setdefault(smiles, [])
                if chebi_id not in ids:
                    ids.append(chebi_id)
            if inchikey_value:
                inchikey_key = inchikey_value.upper()
                ids = inchikey_candidates.setdefault(inchikey_key, [])
                if chebi_id not in ids:
                    ids.append(chebi_id)

    smiles_to_chebi, smiles_collisions = _resolve_lookup_collisions(
        smiles_candidates,
        "SMILES",
    )
    inchikey_to_chebi, inchikey_collisions = _resolve_lookup_collisions(
        inchikey_candidates,
        "InChIKey",
    )

    return smiles_to_chebi, inchikey_to_chebi, smiles_collisions, inchikey_collisions


_local_maps = None


def _get_local_maps():
    """Lazily load and memoize the local SMILES/InChIKey -> ChEBI lookup tables.

    Loaded on first use rather than at import time, so importing this module doesn't
    require the data folder to be configured yet (see chebin.config.set_data_dir).
    """
    global _local_maps
    if _local_maps is None:
        _local_maps = _load_local_smiles_and_inchikey_maps()
    return _local_maps


def convert_smiles_to_chebi(smiles_string, use_parents=False):
    """Convert a single SMILES string to ChEBI IDs.

    Returns (chebi_ids_list, was_resolved, ambiguous_match). ambiguous_match is
    None unless the matched SMILES/InChIKey is asserted by more than one ChEBI
    term in the local lookup table, in which case it's (chosen_chebi_id,
    all_chebi_ids) so the caller can surface the ambiguity to the user.

    use_parents: if no direct ChEBI ID can be found (local table or remote
    lookup), fall back to a remote classification call and use its direct
    parent ChEBI IDs instead of leaving the SMILES unresolved.
    """
    (
        local_smiles_to_chebi,
        local_inchikey_to_chebi,
        local_smiles_collisions,
        local_inchikey_collisions,
    ) = _get_local_maps()

    chebi_ids = []
    was_resolved = False
    ambiguous_match = None
    cleaned_smiles = _clean_lookup_value(smiles_string)
    try:
        mol = Chem.MolFromSmiles(cleaned_smiles)
    except Exception as error:  # noqa: BLE001
        print(f"Warning: failed to parse SMILES {cleaned_smiles}: {error}")
        mol = None
    canonical_smiles = Chem.MolToSmiles(mol) if mol is not None else None

    # First try a direct SMILES lookup against the local leaf-class table, comparing
    # canonical forms so the match doesn't depend on how either SMILES was written.
    lookup_key = canonical_smiles or cleaned_smiles
    local_chebi_id = local_smiles_to_chebi.get(lookup_key)
    if local_chebi_id:
        chosen_id = local_chebi_id.replace("_", ":", 1)
        chebi_ids.append(chosen_id)
        was_resolved = True
        if lookup_key in local_smiles_collisions:
            ambiguous_match = (
                chosen_id,
                [
                    candidate.replace("_", ":", 1)
                    for candidate in local_smiles_collisions[lookup_key]
                ],
            )
        print(
            f"Found local ChEBI ID from exact SMILES match: {chosen_id} for SMILES {cleaned_smiles}",
        )
        return chebi_ids, was_resolved, ambiguous_match

    # If no direct SMILES match exists, try InChIKey -> ChEBI using RDKit to compute
    # the InChIKey for the submitted SMILES.
    try:
        if mol is not None:
            user_inchikey = inchi.MolToInchiKey(mol).upper()
            local_chebi_id = local_inchikey_to_chebi.get(user_inchikey)
            if local_chebi_id:
                chosen_id = local_chebi_id.replace("_", ":", 1)
                chebi_ids.append(chosen_id)
                was_resolved = True
                if user_inchikey in local_inchikey_collisions:
                    ambiguous_match = (
                        chosen_id,
                        [
                            candidate.replace("_", ":", 1)
                            for candidate in local_inchikey_collisions[user_inchikey]
                        ],
                    )
                print(
                    f"Found local ChEBI ID from InChIKey match: {chosen_id} "
                    f"for SMILES {cleaned_smiles} (InChIKey {user_inchikey})",
                )
                return chebi_ids, was_resolved, ambiguous_match
    except Exception as error:  # noqa: BLE001
        print(
            f"Warning: failed to compute InChIKey for SMILES {cleaned_smiles}: {error}",
        )

    # Get details from ChEBI lookup to check for a direct match to a ChEBI ID
    response = requests.post(
        CHEBIFIER_DETAILS_URL,
        json={
            "type": "type",
            "smiles": cleaned_smiles,
            "selectedModels": {
                "ChEBI Lookup": True,
            },
        },
    )

    lookup_model = response.json().get("models", {}).get("ChEBI Lookup", {})

    # Prefer the API's structured chebi_ids list, falling back to pulling the IDs
    # out of the human-readable highlights text (as of 2026-08 the API's chebi_ids
    # field isn't reliably populated, so this fallback is load-bearing).
    lookup_ids = lookup_model.get("chebi_ids")
    if not lookup_ids:
        lookup_infotext = lookup_model.get("highlights", [])
        if lookup_infotext:
            lookup_ids = re.findall(r"CHEBI:(\d+)", lookup_infotext[0][1])

    if lookup_ids:
        # As with the local lookups above, an ambiguous hit contributes a single
        # node and reports the runners-up as alternatives rather than adding every
        # candidate to the study set. The lowest ChEBI ID is used as the tie-break
        # so that every ambiguous match resolves the same way, whichever table or
        # lookup produced it.
        matched_ids = sorted(
            (f"CHEBI:{chebi_id}" for chebi_id in lookup_ids),
            key=lambda cid: int(cid.split(":")[1]),
        )
        chebi_ids.append(matched_ids[0])
        was_resolved = True
        if len(matched_ids) > 1:
            ambiguous_match = (matched_ids[0], matched_ids[1:])
        print(
            f"Found ChEBI ID from lookup: {matched_ids[0]} for SMILES {cleaned_smiles}"
            + (
                f" (ambiguous, also matched {', '.join(matched_ids[1:])})"
                if len(matched_ids) > 1
                else ""
            ),
        )
    elif use_parents:
        print(
            f"No direct ChEBI ID found from lookup for SMILES {cleaned_smiles}, attempting classification...",
        )
        response = requests.post(
            CHEBIFIER_CLASSIFY_URL,
            json={
                "smiles": cleaned_smiles,
                "ontology": False,
                "selectedModels": {
                    "ELECTRA (ChEBI50-3STAR)": True,
                },
            },
        )

        direct_parents = response.json().get("direct_parents")
        if direct_parents:
            # Extract ChEBI IDs from all parent lists
            for parent_list in direct_parents:
                if parent_list is not None:
                    parent_ids = [f"CHEBI:{parent[0]}" for parent in parent_list]
                    chebi_ids.extend(parent_ids)
                    print(
                        f"Found direct parent ChEBI IDs from classification for SMILES {cleaned_smiles}: {parent_ids}",
                    )
                else:
                    print(
                        f"No parents found in one of the classification results for SMILES {cleaned_smiles}",
                    )
                    print(
                        f"Classification response content: {preview(response.content)}",
                    )
            if chebi_ids:
                was_resolved = True

            print(
                f"Found {len(chebi_ids)} ChEBI IDs from classification for SMILES {cleaned_smiles}",
            )

    else:
        print(
            f"No direct ChEBI ID found from lookup for SMILES {cleaned_smiles}, excluding from analysis.",
        )

    return chebi_ids, was_resolved, ambiguous_match


def is_smiles(value: str) -> bool:
    """Heuristic: does this string look like a SMILES rather than a ChEBI ID?

    Not a SMILES if it's already a ChEBI ID in any of the forms
    :mod:`chebin.calculations.chebi_ids` recognises (``CHEBI:12345``,
    ``chebi:12345``, ``CHEBI_12345``, the bare number, ...) or an IRI; otherwise a
    SMILES typically contains lowercase letters, parentheses, or other structural
    characters a bare ChEBI ID never does. The ChEBI check has to come first: the
    letters in ``chebi`` are themselves SMILES atoms, so a lowercased prefix would
    otherwise be read as a structure.
    """
    value = value.strip()
    if looks_like_chebi_id(value) or value.startswith(("http://", "https://")):
        return False
    smiles_chars = set("cCnNoOpPsSFfIiBbr[]()=#@+-\\/")
    return any(c in smiles_chars for c in value)


def _record_ambiguous(ambiguous_matches, smiles, ambiguous_match):
    if ambiguous_match is None:
        return
    chosen, all_ids = ambiguous_match
    ambiguous_matches.append(
        {
            "smiles": smiles,
            "chosen": chosen,
            "alternatives": [cid for cid in all_ids if cid != chosen],
        },
    )


def smiles_list_to_studyset(smiles_list, use_parents=False):
    """Resolve a list of SMILES strings (ChEBI IDs allowed too) to a ChEBI-ID study set.

    Entries that are already ChEBI IDs (``CHEBI:12345``/``CHEBI_12345`` or a full
    IRI, per is_smiles()) are passed through unchanged rather than run through the
    SMILES lookup -- so a list can freely mix SMILES and ChEBI IDs.

    Returns (studyset_list, unresolved_smiles, ambiguous_matches). studyset_list is
    ready to pass to run_enrichment_analysis and friends -- they normalize ChEBI ID
    strings to full IRIs themselves via normalize_id().
    """
    studyset_list = []
    unresolved_smiles = []
    ambiguous_matches = []

    for entry in smiles_list:
        if not is_smiles(entry):
            studyset_list.append(entry)
            continue
        chebi_ids, was_resolved, ambiguous_match = convert_smiles_to_chebi(
            entry,
            use_parents=use_parents,
        )
        if not was_resolved:
            unresolved_smiles.append(entry)
        _record_ambiguous(ambiguous_matches, entry, ambiguous_match)
        studyset_list.extend(chebi_ids)

    return studyset_list, unresolved_smiles, ambiguous_matches


def smiles_weights_to_chebi_weights(smiles_weights, use_parents=False):
    """Resolve a {SMILES: weight} dict (ChEBI IDs allowed too) to a {ChEBI ID: weight} dict.

    Keys that are already ChEBI IDs (per is_smiles()) are passed through unchanged
    rather than run through the SMILES lookup -- so a dict can freely mix SMILES
    and ChEBI ID keys. Each SMILES's weight is applied to every ChEBI ID it
    resolves to. Returns (weights_dict, unresolved_smiles, ambiguous_matches).
    """
    weights_dict = {}
    unresolved_smiles = []
    ambiguous_matches = []

    for entry, weight in smiles_weights.items():
        if not is_smiles(entry):
            weights_dict[entry] = weight
            continue
        chebi_ids, was_resolved, ambiguous_match = convert_smiles_to_chebi(
            entry,
            use_parents=use_parents,
        )
        if not was_resolved:
            unresolved_smiles.append(entry)
        _record_ambiguous(ambiguous_matches, entry, ambiguous_match)
        for chebi_id in chebi_ids:
            weights_dict[chebi_id] = weight

    return weights_dict, unresolved_smiles, ambiguous_matches
