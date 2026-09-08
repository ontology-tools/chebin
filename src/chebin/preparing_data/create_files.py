"""Create all derived data files for the project."""

import json
import os
import shutil
import time
from datetime import UTC
from pathlib import Path

from rdkit import RDLogger  # type: ignore

from chebin.calculations.pre_fishers_calculations import build_class_to_leaf_map
from chebin.calculations.prepare_role_calculations import (
    create_class_to_all_roles_map,
    create_leaves_to_all_roles_map,
    create_roles_to_all_leaves_map,
    find_has_role_connections_from_owl,
)
from chebin.preparing_data.BiGG.get_model import (
    download_model_json,
    gather_recon3d_leaves,
)
from chebin.preparing_data.hmdb.extract_hmdb import extract_hmdb_to_file
from chebin.preparing_data.hmdb.filter_hmdb_statuses import (
    filter_hmdb_extract_by_status,
)
from chebin.preparing_data.load_chebi import load_chebi
from chebin.preparing_data.pruning_smiles import (
    build_parent_map,
    find_leaf_classes_with_smiles_and_deprecated,
    map_names_to_classes,
    save_leaf_classes_with_smiles,
)
from chebin.preparing_data.pruning_split_up_structure import (
    identify_structural_vs_functional,
)
from chebin.preparing_data.wikidata.combine_human_datasets import combine_datasets
from chebin.preparing_data.wikidata.find_missing_chebis import run_find_missing_chebis
from chebin.preparing_data.wikidata.get_inchikeys import convert_smiles_file
from chebin.preparing_data.wikidata.get_lotus import (
    download_lotus_arabidopsis_thaliana,
    download_lotus_homo_sapiens,
)
from chebin.preparing_data.wikidata.get_wikidata_lotus import (
    connect_lotus_csv_to_chebi_ids,
)
from chebin.preparing_data.wikidata.narrow_background_fishers import (
    gather_narrow_leaves,
)


def _run_stage(stage_name, stage_timings, func, *args, **kwargs):
    """Run a stage, record elapsed time, and print a compact timing line."""
    stage_start = time.time()
    result = func(*args, **kwargs)
    elapsed = time.time() - stage_start
    stage_timings.append((stage_name, elapsed))
    print(f"[TIMING] {stage_name}: {elapsed:.2f}s ({elapsed / 60:.2f} min)")
    return result


def _print_timing_summary(stage_timings, total_elapsed):
    """Print stage timings and percentages of the total runtime."""
    print("\n" + "=" * 80)
    print("PIPELINE TIMING SUMMARY")
    print("=" * 80)

    if total_elapsed <= 0:
        total_elapsed = 1e-9

    sorted_timings = sorted(stage_timings, key=lambda x: x[1], reverse=True)
    for stage_name, elapsed in sorted_timings:
        pct = (elapsed / total_elapsed) * 100
        print(
            f"{stage_name}: {elapsed:.2f}s ({elapsed / 60:.2f} min, {pct:.2f}% of total)",
        )

    measured_total = sum(elapsed for _, elapsed in stage_timings)
    overhead = total_elapsed - measured_total
    overhead_pct = (overhead / total_elapsed) * 100
    print("-" * 80)
    print(
        f"Untracked/overhead: {overhead:.2f}s ({overhead / 60:.2f} min, {overhead_pct:.2f}% of total)",
    )
    print("=" * 80)


#: Subfolder for inputs that were downloaded or supplied by the user.
SOURCE_SUBDIR = "source_files"
#: Subfolder for files produced along the way that no calculation ever reads.
INTERMEDIATE_SUBDIR = "intermediate_files"


def find_input_file(data_folder, filename):
    """Locate an input file in the data folder or any of its immediate subfolders.

    The HMDB XML is the one input that cannot be fetched automatically -- the user
    downloads it by hand -- so it is worth finding wherever they reasonably put it,
    whether that is the data folder itself or `source_files/` inside it.

    Returns:
        str | None: Path to the file, or None if it isn't anywhere we looked.
    """
    direct = os.path.join(data_folder, filename)
    if os.path.exists(direct):
        return direct

    if os.path.isdir(data_folder):
        for entry in sorted(os.scandir(data_folder), key=lambda e: e.name):
            if entry.is_dir():
                nested = os.path.join(entry.path, filename)
                if os.path.exists(nested):
                    return nested

    return None


SOURCE_README = """\
# Source files -- not needed for calculations

Inputs the pipeline downloaded or that you supplied. **No enrichment analysis reads
anything in this folder**, so you can delete it to reclaim space.

What deleting each file costs:

- `chebi.owl`, `Recon3D.json`, `lotus_*.csv` -- nothing. Every run of `create_all_files()`
  re-downloads them regardless, so keeping them here saves no work. (ChEBI in particular
  is deliberately re-fetched each time, because it is updated continuously.)
- `hmdb_metabolites.xml` -- **this is the one to keep.** It cannot be downloaded
  automatically. If you delete it you have to fetch 'All Metabolites' again by hand from
  https://hmdb.ca/downloads, and until you do, the first Homo sapiens background
  (`human_entities_leaves.json`) is skipped on the next run.

The pipeline looks for `hmdb_metabolites.xml` here and in the data folder itself, so
either location works.
"""

INTERMEDIATE_README = """\
# Intermediate files -- not needed for calculations

Working files produced while building the data folder: each one is written by a stage,
consumed by the next, and then never read again. **No enrichment analysis reads anything
in this folder**, so it is safe to delete.

Deleting it costs nothing unless you want the files back, which means re-running
`create_all_files()`.
"""


def write_folder_readme(folder, text):
    """Drop a README into a generated subfolder explaining what it is."""
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "README.md"), "w", encoding="utf-8") as f:
        f.write(text)


def rename_folder(old_name, new_name):
    """
    Rename a folder from old_name to new_name.

    Args:
        old_name (str): Current folder name
        new_name (str): New folder name

    Returns:
        bool: True if rename was successful, False otherwise
    """
    try:
        if os.path.exists(old_name):
            Path(old_name).rename(new_name)
            print(f"Folder renamed: '{old_name}' -> '{new_name}'")
            return True
        else:
            print(f"Folder '{old_name}' does not exist")
            return False
    except OSError as e:
        print(f"Error renaming folder: {e}")
        return False


def create_temp_data_folder(new_data_folder="data_new"):
    """
    Create a temporary folder for new files.
    """
    if os.path.exists(f"{new_data_folder}"):
        print(f"Warning: '{new_data_folder}' already exists. Removing it.")
        shutil.rmtree(f"{new_data_folder}")
    os.makedirs(f"{new_data_folder}")
    print(f"Created temporary '{new_data_folder}' folder")


def finalize_folder_structure(new_data_folder="data_new", old_data_folder="data"):
    """
    After all files are created, rename folders:
    1. data -> data_last_used_YYYY.MM.DD (or with counter if needed)
    2. data_new -> data
    """
    from datetime import datetime

    timestamp = datetime.now(UTC).strftime("%Y.%m.%d")

    # Rename old data folder if it exists
    if os.path.exists(f"{old_data_folder}"):
        base_name = f"{old_data_folder}_last_used_{timestamp}"
        old_data_name = base_name
        counter = 1
        while os.path.exists(old_data_name):
            old_data_name = f"{base_name}_{counter}"
            counter += 1
        rename_folder(f"{old_data_folder}", old_data_name)

    # Rename the temporary folder into place
    if os.path.exists(f"{new_data_folder}"):
        rename_folder(f"{new_data_folder}", f"{old_data_folder}")


def cleanup_old_data_folders(old_data_folder="data", max_folders=3):
    """
    Keep only the most recent 'max_folders' of old data folders and delete the rest.

    Args:
        old_data_folder (str): Name of the old data folder
        max_folders (int): Maximum number of old data folders to keep
    """
    import re

    # The backups are written as siblings of the data folder itself, so look for them
    # there rather than in the current working directory.
    data_path = Path(old_data_folder)
    parent_dir = data_path.parent
    base_name = data_path.name

    # Match folders named "data_last_used_YYYY.MM.DD" or "data_last_used_YYYY.MM.DD_N"
    pattern = re.compile(
        rf"{re.escape(base_name)}_last_used_(\d{{4}}\.\d{{2}}\.\d{{2}})(?:_(\d+))?$",
    )

    old_data_folders = []
    for entry in parent_dir.iterdir():
        if not entry.is_dir():
            continue
        match = pattern.match(entry.name)
        if match:
            date_str, counter = match.groups()
            old_data_folders.append((entry, date_str, int(counter) if counter else 0))

    # Sort by the date (and counter, for same-day folders) encoded in the name, newest first
    old_data_folders.sort(key=lambda item: (item[1], item[2]), reverse=True)

    # Keep only the most recent 'max_folders'
    folders_to_delete = old_data_folders[max_folders:]

    for folder, _, _ in folders_to_delete:
        print(f"Deleting old data folder: {folder}")
        shutil.rmtree(folder)


def create_all_files_with_backup(data_folder="data"):
    """
    Create all derived data files for the project, with a backup of the old data folder. Observe that this function can take up to an hour to run.
    """
    data_folder_new = f"{data_folder}_new"

    ### Create temporary data folder for new files
    print("Creating temporary data folder...")
    create_temp_data_folder(new_data_folder=data_folder_new)

    # Carry the HMDB XML over to the new folder before the rename, if we have it.
    # It is optional: without it, create_all_files simply skips the first human
    # background (see the warning it prints).
    # The HMDB XML is the only input worth carrying across: ChEBI, Recon3D and the LOTUS
    # CSVs are re-downloaded on every run anyway, but this one was placed by hand and
    # cannot be fetched automatically. It may sit in the old folder itself or in a
    # subfolder, so look in both.
    hmdb_xml_src = find_input_file(data_folder, "hmdb_metabolites.xml")
    if hmdb_xml_src:
        hmdb_xml_dst = f"{data_folder_new}/{SOURCE_SUBDIR}/hmdb_metabolites.xml"
        os.makedirs(os.path.dirname(hmdb_xml_dst), exist_ok=True)
        print(f"Copying HMDB XML to '{data_folder_new}/{SOURCE_SUBDIR}'...")
        shutil.copy2(hmdb_xml_src, hmdb_xml_dst)
        print(f"Copied {hmdb_xml_src} to {hmdb_xml_dst}")
    else:
        print(
            f"hmdb_metabolites.xml not found in '{data_folder}' or its subfolders; "
            f"skipping the HMDB copy.",
        )

    create_all_files(data_folder=data_folder_new)

    ### Finalize folder structure: rename old data and move data_new to data
    print("Finalizing folder structure...")
    finalize_folder_structure(
        new_data_folder=data_folder_new,
        old_data_folder=data_folder,
    )

    ### Clean up old data folders, keeping only the most recent ones
    print("Cleaning up old data folders...")
    cleanup_old_data_folders(old_data_folder=data_folder, max_folders=3)


def create_all_files(data_folder="data"):
    """
    Create all derived data files for the project. Observe that this function can take up to an hour to run.
    Args:
        data_folder (str): Path to the folder where the data files will be created. Defaults to "data".
    """

    # ChEBI's SMILES routinely trip RDKit's "Omitted undefined stereo", "Can't
    # kekulize mol" and "Proton(s) added/removed" warnings, which the pipeline
    # already handles by falling back or skipping. Left on, they produced ~11 MB
    # of stderr per run -- 30x the size of the actual log. Scoped to this
    # function so importing chebin elsewhere (the website) keeps the warnings.
    RDLogger.DisableLog("rdApp.*")

    stage_timings = []
    start_time = time.time()

    ### Two subfolders keep the data folder itself down to just the files the
    ### calculations read; everything else is sorted into one of these, each carrying a
    ### README saying it is safe to delete. Created up front because not every stage
    ### makes its own output directory (prepare_role_calculations does not).
    source_dir = f"{data_folder}/{SOURCE_SUBDIR}"
    intermediate_dir = f"{data_folder}/{INTERMEDIATE_SUBDIR}"
    os.makedirs(source_dir, exist_ok=True)
    os.makedirs(intermediate_dir, exist_ok=True)

    # The HMDB XML is placed by hand, so accept it in the data folder or any subfolder.
    hmdb_xml_file = find_input_file(data_folder, "hmdb_metabolites.xml")
    hmdb_xml_exists = hmdb_xml_file is not None

    # Warn the user already if HMDB XML is missing, so they can download it before the long process starts
    if not hmdb_xml_exists:
        print(
            f"hmdb_metabolites.xml not found in '{data_folder}' or its subfolders. "
            f"The process will continue unless you stop it. "
            f"This file is only needed for one of the human backgrounds. "
            f"To be able to use it, please download it from https://hmdb.ca/downloads "
            f"and place it in '{source_dir}'.",
        )
    else:
        print(f"Using HMDB XML at {hmdb_xml_file}")
        # Found, but not where the rest of the inputs live. Suggest rather than move:
        # it is the user's file and it is large, so relocating it silently would be rude.
        if os.path.dirname(os.path.abspath(hmdb_xml_file)) == os.path.abspath(
            data_folder,
        ):
            print(
                f"  (tip: moving it into '{source_dir}' keeps the data folder itself to "
                f"just the files the calculations need -- both locations work)",
            )

    ### Define properties and file paths
    smiles_property = "https://w3id.org/chemrof/smiles_string"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"
    has_role_property = "http://purl.obolibrary.org/obo/RO_0000087"

    # --- source files: downloaded or user-supplied, read only during the build
    chebi_file = f"{source_dir}/chebi.owl"
    recon3d_json = f"{source_dir}/Recon3D.json"
    lotus_hs_csv = f"{source_dir}/lotus_homo_sapiens.csv"
    lotus_at_csv = f"{source_dir}/lotus_arabidopsis_thaliana.csv"

    # --- intermediates: written by one stage, consumed by the next, then dead
    subclass_map_file = f"{intermediate_dir}/chebi_subclass_map.json"
    roles_map_json = f"{intermediate_dir}/class_to_direct_roles_map.json"  # output file for roles map
    leaves_to_all_roles_json = f"{intermediate_dir}/removed_leaf_classes_to_ALL_roles_map.json"  # output file for leaves to all roles map

    lotus_hs_chebi_tsv = f"{intermediate_dir}/lotus_homo_sapiens_with_chebi_ids.tsv"
    lotus_at_chebi_tsv = (
        f"{intermediate_dir}/lotus_arabidopsis_thaliana_with_chebi_ids.tsv"
    )
    lotus_hs_updated_tsv = (
        f"{intermediate_dir}/lotus_homo_sapiens_with_chebi_ids_updatedchebis.tsv"
    )
    lotus_at_updated_tsv = f"{intermediate_dir}/lotus_arabidopsis_thaliana_with_chebi_ids_updatedchebis.tsv"

    hmdb_extract_tsv = f"{intermediate_dir}/hmdb_metabolites_extract.tsv"
    hmdb_filtered_tsv = (
        f"{intermediate_dir}/hmdb_metabolites_extract_quantified_detected.tsv"
    )
    hmdb_updated_tsv = f"{intermediate_dir}/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv"
    combined_human_tsv = f"{intermediate_dir}/combined_hmdb_wikidata.tsv"

    # --- runtime files: everything below stays at the top level of the data folder,
    # --- because these are what the enrichment analyses and the website actually read.
    leaf_parents_map_file = (
        f"{data_folder}/removed_leaf_classes_to_ALL_parents_map.json"
    )
    removed_leaf_classes_file = f"{data_folder}/removed_leaf_classes_with_smiles.csv"
    leaves_to_all_parents_json = (
        f"{data_folder}/removed_leaf_classes_to_ALL_parents_map.json"
    )
    parent_map_json = f"{data_folder}/chebi_parent_map.json"
    roles_to_all_leaves_json = f"{data_folder}/roles_to_leaves_map.json"  # output file for roles to all leaves map
    class_to_all_roles_json = f"{data_folder}/class_to_all_roles_map.json"  # output file for class to all roles map
    id_to_name_map_json = f"{data_folder}/chebi_id_to_name_map.json"

    leaf_to_ancestors_file = (
        f"{data_folder}/removed_leaf_classes_to_ALL_parents_map.json"
    )
    class_to_leaf_output_file = f"{data_folder}/class_to_leaf_descendants_map.json"

    # Lookup tables used to match compounds to ChEBI IDs. These must be passed
    # explicitly to every stage below: the underlying functions default to "data/...",
    # which would silently read the previous run's files when data_folder is a
    # temporary folder (see create_all_files_with_backup).
    leaves_smiles_csv = removed_leaf_classes_file
    leaves_inchikeys_csv = f"{data_folder}/removed_leaf_classes_with_inchikeys.csv"

    ### Download and load the ChEBI ontology
    print("Downloading and loading ChEBI ontology...")
    # Option: download directly to data_new/ to keep separate from old versions
    chebi_ontology = _run_stage(
        "download_and_load_chebi",
        stage_timings,
        load_chebi,
        download_dir=source_dir,
    )

    ### Find and filter leaf classes with SMILES and deprecated classes, and save the filtered ontology
    print("Removing leaf classes with SMILES and deprecated classes...")
    classes_with_smiles, deprecated_classes = _run_stage(
        "find_leaf_classes_with_smiles_and_deprecated",
        stage_timings,
        find_leaf_classes_with_smiles_and_deprecated,
        chebi_ontology,
        smiles_property,
        deprecated_property,
        subclass_map_file,
        leaf_parents_map_file,
        use_found_leaf_classes=False,
        removed_leaf_classes_file=removed_leaf_classes_file,
    )

    ### Build map of all classes to their direct parents and save as JSON
    print("Building parent map...")
    _run_stage(
        "build_parent_map",
        stage_timings,
        build_parent_map,
        chebi_ontology,
        parent_map_json,
        deprecated_property,
        subclass_map_file=subclass_map_file,
        precomputed_deprecated_classes=deprecated_classes,
    )

    ### Build map with CHEBI short IDs to names and save as JSON
    print("Building ID to name map...")
    _run_stage(
        "map_names_to_classes",
        stage_timings,
        map_names_to_classes,
        chebi_file,
        id_to_name_map_json,
    )

    ### Map all classes to their direct "has_role" connections and save as JSON,
    ### then create a map of all leaf classes to all their "has_role" connections (not just direct ones) and save as JSON
    ### and the reverse map of all roles to all leaf classes that have that role (directly or indirectly) and save as JSON
    print("Building roles maps...")
    _run_stage(
        "find_has_role_connections_from_owl",
        stage_timings,
        find_has_role_connections_from_owl,
        chebi_file,
        has_role_property,
        deprecated_property,
        roles_map_json,
    )
    _run_stage(
        "create_leaves_to_all_roles_map",
        stage_timings,
        create_leaves_to_all_roles_map,
        roles_map_json,
        leaves_to_all_parents_json,
        leaves_to_all_roles_json,
        parent_map_json,
    )
    _run_stage(
        "create_roles_to_all_leaves_map",
        stage_timings,
        create_roles_to_all_leaves_map,
        leaves_to_all_roles_json,
        roles_to_all_leaves_json,
    )

    ### Map each class to its direct roles, roles inherited from ancestor classes, and descendants of those roles in the role hierarchy.
    _run_stage(
        "create_class_to_all_roles_map",
        stage_timings,
        create_class_to_all_roles_map,
        roles_map_json,
        parent_map_json,
        class_to_all_roles_json,
    )

    ### Identify structural vs functional classes (fast path: traverse the precomputed
    ### subclass map instead of the slow per-node ontology API). The flattened map written by
    ### find_leaf_classes gives identical descendant sets for the (non-leaf) roots.
    print("Identifying structural vs functional classes...")
    with open(subclass_map_file) as f:
        subclass_map_for_classification = json.load(f)
    structural_classes, functional_classes, unknown_classes = _run_stage(
        "identify_structural_vs_functional",
        stage_timings,
        identify_structural_vs_functional,
        chebi_ontology,
        subclass_map_for_classification,
    )

    print(f"Identified {len(structural_classes)} structural classes")
    print(f"Identified {len(functional_classes)} functional classes")
    print(f"Identified {len(unknown_classes)} unknown classes")

    ### Save a CSV file with the removed leaf classes (using in-memory class sets)
    print("Saving CSV with removed leaf classes...")
    _run_stage(
        "save_leaf_classes_with_smiles",
        stage_timings,
        save_leaf_classes_with_smiles,
        classes_with_smiles,
        chebi_ontology,
        smiles_property,
        removed_leaf_classes_file,
        structural_classes,
        functional_classes,
        owl_file=chebi_file,  # OPTIMIZATION: Pass OWL file for fast XML parsing of SMILES
        parent_map_file=parent_map_json,  # OPTIMIZATION: Use already-created parent map instead of ontology API
    )

    ### Build a JSON map from each class IRI to ALL its leaf descendants
    print("Building class to leaf descendants map...")
    _run_stage(
        "build_class_to_leaf_map",
        stage_timings,
        build_class_to_leaf_map,
        leaf_to_ancestors_file,
        class_to_leaf_output_file,
    )

    ### Convert SMILES to InChIKeys for the removed leaf classes (needed by find_missing_chebis and the website)
    print("Generating InChIKeys for removed leaf classes...")
    _run_stage(
        "convert_smiles_to_inchikeys",
        stage_timings,
        convert_smiles_file,
        leaves_smiles_csv,
        leaves_inchikeys_csv,
    )

    ### The first human background combines HMDB with LOTUS Homo sapiens, so without the
    ### HMDB XML the whole LOTUS Homo sapiens chain would produce files nothing reads.
    ### Skip it all together rather than paying for the download and the matching.
    if not hmdb_xml_exists:
        print(
            "HMDB XML file not found. Skipping the first human background "
            "(HMDB + LOTUS Homo sapiens) and everything feeding into it.",
        )
    else:
        print("Downloading LOTUS Homo sapiens explorer CSV...")
        _run_stage(
            "download_lotus_homo_sapiens",
            stage_timings,
            download_lotus_homo_sapiens,
            lotus_hs_csv,
        )

        print("Matching ChEBI IDs for LOTUS Homo sapiens compounds...")
        _run_stage(
            "connect_lotus_homo_sapiens_to_chebi_ids",
            stage_timings,
            connect_lotus_csv_to_chebi_ids,
            lotus_hs_csv,
            lotus_hs_chebi_tsv,
            inchikeys_csv=leaves_inchikeys_csv,
            smiles_csv=leaves_smiles_csv,
        )

        print("Creating HMDB output file...")
        _run_stage(
            "create_hmdb_output_file",
            stage_timings,
            extract_hmdb_to_file,
            hmdb_xml_file,
            hmdb_extract_tsv,
        )

        print("Filtering HMDB statuses...")
        _run_stage(
            "filter_hmdb_statuses",
            stage_timings,
            filter_hmdb_extract_by_status,
            input_file=hmdb_extract_tsv,
            output_file=hmdb_filtered_tsv,
        )

        print("Finding missing ChEBI IDs for LOTUS (Homo sapiens)...")
        _run_stage(
            "find_missing_chebis_lotus_hs",
            stage_timings,
            run_find_missing_chebis,
            "lotus_hs",
            compounds_file=lotus_hs_chebi_tsv,
            output_file=lotus_hs_updated_tsv,
        )

        print("Finding missing ChEBI IDs for HMDB...")
        _run_stage(
            "find_missing_chebis_hmdb",
            stage_timings,
            run_find_missing_chebis,
            "hmdb",
            compounds_file=hmdb_filtered_tsv,
            output_file=hmdb_updated_tsv,
        )

        print("Combining human datasets...")
        _run_stage(
            "combine_datasets",
            stage_timings,
            combine_datasets,
            hmdb_path=hmdb_updated_tsv,
            wikidata_path=lotus_hs_updated_tsv,
            output_path=combined_human_tsv,
        )

        print("Gathering narrow leaf classes (Homo sapiens)...")
        _run_stage(
            "gather_narrow_leaves_homo_sapiens",
            stage_timings,
            gather_narrow_leaves,
            compounds_tsv=combined_human_tsv,
            leaves_csv=leaves_smiles_csv,
            class_to_leaf_map=class_to_leaf_output_file,
            output_json=f"{data_folder}/human_entities_leaves.json",
            taxon_label="homo_sapiens",
        )

    print("Downloading LOTUS Arabidopsis thaliana explorer CSV...")
    _run_stage(
        "download_lotus_arabidopsis_thaliana",
        stage_timings,
        download_lotus_arabidopsis_thaliana,
        lotus_at_csv,
    )

    print("Matching ChEBI IDs for LOTUS Arabidopsis thaliana compounds...")
    _run_stage(
        "connect_lotus_arabidopsis_thaliana_to_chebi_ids",
        stage_timings,
        connect_lotus_csv_to_chebi_ids,
        lotus_at_csv,
        lotus_at_chebi_tsv,
        inchikeys_csv=leaves_inchikeys_csv,
        smiles_csv=leaves_smiles_csv,
    )

    print("Finding missing ChEBI IDs for LOTUS (Arabidopsis thaliana)...")
    _run_stage(
        "find_missing_chebis_lotus_at",
        stage_timings,
        run_find_missing_chebis,
        "lotus_at",
        compounds_file=lotus_at_chebi_tsv,
        output_file=lotus_at_updated_tsv,
    )

    print("Gathering narrow leaf classes (Arabidopsis thaliana)...")
    _run_stage(
        "gather_narrow_leaves_arabidopsis_thaliana",
        stage_timings,
        gather_narrow_leaves,
        compounds_tsv=lotus_at_updated_tsv,
        leaves_csv=leaves_smiles_csv,
        class_to_leaf_map=class_to_leaf_output_file,
        output_json=f"{data_folder}/arabidopsis_thaliana_leaves.json",
        taxon_label="arabidopsis_thaliana",
    )

    print("Downloading Recon3D model...")
    _run_stage(
        "download_recon3d_model",
        stage_timings,
        download_model_json,
        "Recon3D",
        recon3d_json,
    )

    print("Gathering narrow leaf classes (Endogenous human / Recon3D)...")
    _run_stage(
        "gather_narrow_leaves_endogenous_human",
        stage_timings,
        gather_recon3d_leaves,
        model_json_path=recon3d_json,
        leaves_csv=leaves_smiles_csv,
        class_to_leaf_map_json=class_to_leaf_output_file,
        output_json=f"{data_folder}/recon3d_leaves.json",
    )

    ### Label the two subfolders so a user opening the data folder can tell at a glance
    ### which files matter and what deleting the rest costs.
    write_folder_readme(source_dir, SOURCE_README)
    write_folder_readme(intermediate_dir, INTERMEDIATE_README)

    print("All files created successfully!")
    print(
        f"  Files needed for calculations: {data_folder}/\n"
        f"  Safe to delete:                {source_dir}/ and {intermediate_dir}/ "
        f"(see the README in each)",
    )
    end_time = time.time()
    elapsed_time = end_time - start_time
    _print_timing_summary(stage_timings, elapsed_time)
    print(
        f"Total execution time: {elapsed_time:.2f} seconds or {elapsed_time / 60:.2f} minutes or {elapsed_time / 3600:.2f} hours",
    )


if __name__ == "__main__":
    create_all_files_with_backup()
