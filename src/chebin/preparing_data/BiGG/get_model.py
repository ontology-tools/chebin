import json
import re
import time
import urllib.error
import urllib.request

import pandas as pd

# ruff: noqa
import requests

from chebin.calculations.fishers_calculations import normalize_id
from chebin.calculations.predicted_parents import LEAF_EXPANSION_LIMIT

MODEL_URL_TEMPLATE = "http://bigg.ucsd.edu/static/models/{model_id}.json"
UNICHEM_URL = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
UNICHEM_HMDB_SOURCE_ID = 18  # UniChem's source ID for HMDB
# Skip expanding a class this generic: the same cutoff as everywhere else.
MAX_LEAF_DESCENDANTS = LEAF_EXPANSION_LIMIT


def normalize_hmdb_id(hmdb_id):
    """Recon3D stores HMDB IDs in the old 5-digit form ('HMDB06537'); UniChem
    expects the modern 7-digit zero-padded form ('HMDB0006537')."""
    digits = re.sub(r"^HMDB", "", hmdb_id.strip())
    if not digits.isdigit():
        return hmdb_id.strip()
    return f"HMDB{int(digits):07d}"


def download_model_json(model_id, out_path):
    url = MODEL_URL_TEMPLATE.format(model_id=model_id)
    resp = requests.get(url)
    resp.raise_for_status()
    with open(out_path, "wb") as f:
        f.write(resp.content)


def load_unique_compounds(model_json_path):
    """Collapse per-compartment metabolite entries (eg '10fthf_c', '10fthf_m') into one record per compound."""
    with open(model_json_path) as f:
        data = json.load(f)

    compounds = {}
    for met in data["metabolites"]:
        base_id = re.sub(r"_[a-z0-9]+$", "", met["id"])
        ann = met.get("annotation", {})
        rec = compounds.setdefault(
            base_id,
            {
                "name": met.get("name", ""),
                "chebi": set(),
                "hmdb": set(),
                "inchi_key": set(),
            },
        )
        rec["chebi"].update(ann.get("chebi", []))
        rec["hmdb"].update(ann.get("hmdb", []))
        rec["inchi_key"].update(ann.get("inchi_key", []))
    return compounds


def _query_unichem(body):
    """POST to UniChem and pull out any ChEBI IDs from the matched structure bucket."""
    req = urllib.request.Request(
        UNICHEM_URL,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            payload = json.load(resp)
    except (urllib.error.URLError, urllib.error.HTTPError, OSError):
        return set()

    chebi_ids = set()
    for compound in payload.get("compounds", []):
        for source in compound.get("sources", []):
            if source.get("shortName") == "chebi":
                chebi_ids.add(source["compoundId"])
    return chebi_ids


def query_unichem_by_hmdb(hmdb_id):
    """Cross-reference an HMDB ID to ChEBI IDs via UniChem (exact structure match)."""
    return _query_unichem(
        {
            "type": "sourceID",
            "compound": normalize_hmdb_id(hmdb_id),
            "sourceID": UNICHEM_HMDB_SOURCE_ID,
        },
    )


def query_unichem_by_inchikey(inchi_key):
    """Cross-reference an InChIKey to ChEBI IDs via UniChem (exact structure match)."""
    return _query_unichem({"type": "inchikey", "compound": inchi_key.strip()})


def resolve_to_leaves(
    candidate_chebi_ids,
    leaf_ids,
    class_to_leaf_map,
    max_descendants=MAX_LEAF_DESCENDANTS,
):
    """Keep candidates that are already leaves; if none are, expand non-leaf
    candidates to their leaf descendants, skipping any class with more than
    `max_descendants` (too generic to usefully narrow the background).

    Returns (resolved_leaves, method) where method is 'leaf', 'expansion', or 'none'.
    """
    candidates = {normalize_id(c) for c in candidate_chebi_ids}

    leaf_hits = candidates & leaf_ids
    if leaf_hits:
        return leaf_hits, "leaf"

    expanded = set()
    for candidate in candidates:
        descendants = class_to_leaf_map.get(candidate, [])
        if descendants and len(descendants) <= max_descendants:
            expanded.update(descendants)
    if expanded:
        return expanded, "expansion"

    return set(), "none"


def gather_recon3d_leaves(
    model_json_path="data/Recon3D.json",
    leaves_csv="data/removed_leaf_classes_with_smiles.csv",
    class_to_leaf_map_json="data/class_to_leaf_descendants_map.json",
    output_json="data/recon3d_leaves.json",
):
    compounds = load_unique_compounds(model_json_path)

    leaf_ids = set(pd.read_csv(leaves_csv)["IRI"].values)
    with open(class_to_leaf_map_json) as f:
        class_to_leaf_map = json.load(f)

    recon3d_leaves = set()
    unresolved = []
    counts = {
        "resolved_leaf_direct": 0,
        "resolved_via_expansion": 0,
        "resolved_via_unichem_inchikey_leaf": 0,
        "resolved_via_unichem_inchikey_expansion": 0,
        "resolved_via_unichem_hmdb_leaf": 0,
        "resolved_via_unichem_hmdb_expansion": 0,
        "unresolved": 0,
    }

    for compound_id, rec in compounds.items():
        candidates = set(rec["chebi"])
        unichem_source = None

        # No ChEBI from BiGG: fall back to UniChem exact-structure lookup.
        # Try InChIKey first (most reliable), then the old-format HMDB IDs.
        if not candidates and rec["inchi_key"]:
            for inchi_key in rec["inchi_key"]:
                candidates.update(query_unichem_by_inchikey(inchi_key))
                time.sleep(0.3)
            if candidates:
                unichem_source = "inchikey"

        if not candidates and rec["hmdb"]:
            for hmdb_id in rec["hmdb"]:
                candidates.update(query_unichem_by_hmdb(hmdb_id))
                time.sleep(0.3)
            if candidates:
                unichem_source = "hmdb"

        if not candidates:
            unresolved.append(compound_id)
            counts["unresolved"] += 1
            continue

        leaves, method = resolve_to_leaves(candidates, leaf_ids, class_to_leaf_map)

        if not leaves:
            unresolved.append(compound_id)
            counts["unresolved"] += 1
            continue

        recon3d_leaves.update(leaves)
        # method is 'leaf' (candidate was itself a leaf) or 'expansion' (candidate
        # was a parent, expanded to its leaf descendants).
        if unichem_source == "inchikey":
            counts[f"resolved_via_unichem_inchikey_{method}"] += 1
        elif unichem_source == "hmdb":
            counts[f"resolved_via_unichem_hmdb_{method}"] += 1
        elif method == "leaf":
            counts["resolved_leaf_direct"] += 1
        else:
            counts["resolved_via_expansion"] += 1

    payload = {
        "taxon_label": "homo_sapiens_recon3d",
        "model_json_path": model_json_path,
        "n_compounds_total": len(compounds),
        "n_compounds_resolved": len(compounds) - len(unresolved),
        "n_compounds_unresolved": len(unresolved),
        "resolution_counts": counts,
        "narrow_leaves": sorted(recon3d_leaves),
    }

    with open(output_json, "w") as f:
        json.dump(payload, f, indent=2)

    print(f"Total compounds (BiGG, compartments collapsed): {len(compounds)}")
    print(
        f"  Resolved directly (BiGG ChEBI ID is already a leaf):       {counts['resolved_leaf_direct']}",
    )
    print(
        f"  Resolved via parent expansion (<= {MAX_LEAF_DESCENDANTS} leaf descendants):   {counts['resolved_via_expansion']}",
    )
    print(
        f"  Resolved via UniChem InChIKey -> already a leaf:           {counts['resolved_via_unichem_inchikey_leaf']}",
    )
    print(
        f"  Resolved via UniChem InChIKey -> expanded to leaves:       {counts['resolved_via_unichem_inchikey_expansion']}",
    )
    print(
        f"  Resolved via UniChem HMDB     -> already a leaf:           {counts['resolved_via_unichem_hmdb_leaf']}",
    )
    print(
        f"  Resolved via UniChem HMDB     -> expanded to leaves:       {counts['resolved_via_unichem_hmdb_expansion']}",
    )
    print(
        f"  Unresolved (left out of the background):                   {counts['unresolved']}",
    )
    print(
        f"Total leaf ChEBI classes in the resulting background: {len(recon3d_leaves)}",
    )
    print(f"Saved to {output_json}.")

    return recon3d_leaves, unresolved


if __name__ == "__main__":
    download_model_json("Recon3D", "data/Recon3D.json")
    gather_recon3d_leaves()
