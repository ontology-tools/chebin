import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import glob
import re
import time
import uuid

from flask import Flask, redirect, render_template, request, session, url_for

from chebin.calculations.chebi_ids import collapse_chebi_prefixes, to_chebi_curie
from chebin.calculations.data_files import NARROW_BACKGROUND_LEAVES
from chebin.calculations.fishers_calculations import (
    run_enrichment_analysis,
    run_enrichment_analysis_plain_enrich_pruning_strategy,
)
from chebin.calculations.smiles_lookup import convert_smiles_to_chebi, is_smiles
from chebin.calculations.visualitations_and_pruning import graph_to_cytospace_json
from chebin.calculations.weighted_calculations import (
    run_weighted_enrichment_analysis,
    run_weighted_enrichment_analysis_plain_enrich_pruning_strategy,
    run_weighted_narrow_background_enrichment_analysis,
    run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy,
)
from chebin.config import data_path, set_data_dir
from chebin.preparing_data.wikidata.narrow_background_fishers import (
    run_narrow_background_enrichment_analysis,
    run_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy,
)

app = Flask(__name__)
app.secret_key = "your_secret_key"  # Replace with a secure secret key

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.argv[0] = os.path.abspath(
    sys.argv[0],
)  # keep Werkzeug's debug-reloader re-exec working after chdir below
os.chdir(BASE_DIR)

# Tell chebin where the data files are, rather than relying on the chdir above:
# the package resolves them against the working directory by default.
set_data_dir(os.path.join(BASE_DIR, "data"))

# Human-readable names for the backgrounds, used when one is unavailable.
BACKGROUND_LABELS = {
    "human": "Homo sapiens-based background 1",
    "arabidopsis_thaliana": "Arabidopsis thaliana-based background",
    "endogenous_human": "Homo sapiens-based background 2",
}


def _missing_background_file(background):
    """Path of the leaves JSON a narrow background needs, if it isn't on disk.

    Some backgrounds are optional in the data pipeline: the first human background is
    only built when the HMDB XML was present, so a local data folder may legitimately
    lack it. Checked per request rather than at import, since the file appears as soon
    as the pipeline is re-run.
    """
    filename = NARROW_BACKGROUND_LEAVES.get(background)
    if filename is None:
        return None
    leaves_json = data_path(filename)
    return None if os.path.exists(leaves_json) else leaves_json


def _render_background_unavailable(background, leaves_json):
    return render_template(
        "background_unavailable.html",
        background_label=BACKGROUND_LABELS.get(background, background),
        leaves_json=leaves_json,
    )


def cleanup_old_graph_files(max_age_hours=24):
    """Remove graph files older than max_age_hours"""
    cutoff = time.time() - (max_age_hours * 3600)
    for filepath in glob.glob("website/static/data/graph_*.json"):
        try:
            if os.path.getmtime(filepath) < cutoff:
                os.remove(filepath)
        except OSError:
            pass


# Load last update time
def get_data_version():
    try:
        with open(os.path.join(BASE_DIR, "data_version.txt")) as f:
            return f.read().strip()
    except FileNotFoundError:
        return "unknown"


@app.context_processor
def inject_data_version():
    return {"data_version": get_data_version()}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/submission", methods=["GET", "POST"])
def submission():
    if request.method == "POST":
        study_set = request.form.get("study_set")
        session["study_set"] = study_set  # Store the study set in session
        classification = request.form.get("classification")
        if classification:
            session["classification"] = (
                classification  # Store classification in session
            )
        smiles_option = request.form.get("smiles_option")
        session["smiles_option"] = smiles_option  # Store smiles option in session
        background = request.form.get("background")
        # Catch an unavailable background here, before the user fills in a study set.
        missing_leaves_json = _missing_background_file(background)
        if missing_leaves_json:
            return _render_background_unavailable(background, missing_leaves_json)
        session["background"] = background  # Store background in session
        # Store user's preference for expanding the human background (checkbox)
        expand_background = bool(request.form.get("expand_background"))
        session["expand_background"] = expand_background

        return render_template("submission.html", user_study_set=study_set)

    return render_template("submission.html", user_study_set=None)


def _line_weight(parts):
    """The weight of a ``<id> <weight>`` study-set line, or None if there isn't one.

    A bare number is a valid ChEBI ID, so a line of two of them (``17079 17080``,
    or ``17079, 17080``) is two IDs rather than an ID and a weight. A second token
    only counts as a weight when the line says so unambiguously: either the first
    token carries a CHEBI prefix, or the weight isn't a whole number.
    """
    if len(parts) < 2:
        return None
    try:
        weight = float(parts[1])
    except ValueError:
        return None
    if parts[0].isdigit() and parts[1].isdigit():
        return None
    return weight


# Used in run_analysis route to parse user input
def parse_studyset(studyset: str):
    # Remove surrounding quotes if present
    studyset = studyset.strip()

    studyset_list = []
    weights_dict = {}
    unresolved_smiles = []
    ambiguous_smiles_matches = []

    if not studyset:
        return studyset_list, weights_dict, unresolved_smiles, ambiguous_smiles_matches

    def normalize_id(raw_id: str) -> str:
        """Fold the ChEBI ID spellings users type into the single CHEBI_12345 form.

        The analysis functions normalize to full IRIs themselves (see
        chebin.calculations.fishers_calculations.normalize_id), so this only has
        to make the ID recognisable to them; anything that isn't a ChEBI ID is
        passed through for them to reject.
        """
        curie = to_chebi_curie(raw_id)
        if curie is not None:
            return curie.replace(":", "_")
        value = raw_id.strip().replace('"', "")
        if value.startswith(("http://", "https://")):
            return value
        return value.replace(":", "_")

    def record_ambiguous(smiles, ambiguous_match):
        if ambiguous_match is None:
            return
        chosen, all_ids = ambiguous_match
        ambiguous_smiles_matches.append(
            {
                "smiles": smiles,
                "chosen": chosen,
                "alternatives": [cid for cid in all_ids if cid != chosen],
            },
        )

    use_parents = session.get("smiles_option") == "use_parents"

    # Split by lines first to support optional weights per line
    for line in studyset.splitlines():
        line = line.strip()
        if not line:
            continue

        # Fold "CHEBI ID: 17079" and friends into "CHEBI:17079" first: the split
        # below would otherwise scatter one ID across several tokens.
        line = collapse_chebi_prefixes(line)

        parts = [p for p in re.split(r"[\s,]+", line) if p]

        weight = _line_weight(parts)
        if weight is not None:
            # Check if first part is SMILES
            if is_smiles(parts[0]):
                chebi_ids, was_resolved, ambiguous_match = convert_smiles_to_chebi(
                    parts[0],
                    use_parents=use_parents,
                )
                if not was_resolved:
                    unresolved_smiles.append(parts[0])
                record_ambiguous(parts[0], ambiguous_match)
                # Apply the same weight to all resulting ChEBI IDs
                for chebi_id in chebi_ids:
                    class_id = normalize_id(chebi_id)
                    studyset_list.append(class_id)
                    weights_dict[class_id] = weight
            else:
                class_id = normalize_id(parts[0])
                studyset_list.append(class_id)
                weights_dict[class_id] = weight
            continue

        # Fallback: treat all parts as IDs without weights
        for part in parts:
            if is_smiles(part):
                chebi_ids, was_resolved, ambiguous_match = convert_smiles_to_chebi(
                    part,
                    use_parents=use_parents,
                )
                if not was_resolved:
                    unresolved_smiles.append(part)
                record_ambiguous(part, ambiguous_match)
                for chebi_id in chebi_ids:
                    class_id = normalize_id(chebi_id)
                    studyset_list.append(class_id)
            else:
                class_id = normalize_id(part)
                studyset_list.append(class_id)

    return studyset_list, weights_dict, unresolved_smiles, ambiguous_smiles_matches


def map_p_value_correction_method(method_name):
    if method_name == "bonferroni":
        return True, False
    elif method_name == "benjamini_hochberg":
        return False, True
    else:
        return False, False  # No correction method selected


@app.route("/run_analysis", methods=["GET", "POST"])
# option to choose pruning methods
def run_analysis():
    # Cleanup old graph files
    cleanup_old_graph_files(max_age_hours=24)

    # Generate or retrieve session ID
    if "session_id" not in session:
        session["session_id"] = str(uuid.uuid4())

    raw_studyset = session.get("study_set")
    if not raw_studyset:
        return redirect(url_for("submission"))
    # Convert multi-line or comma-separated input into a list
    studyset_list, weights_dict, unresolved_smiles, ambiguous_smiles_matches = (
        parse_studyset(raw_studyset)
    )

    # Auto-scale weights if present (only scales up if max < 1000)
    # if weights_dict:
    #     weights_dict = auto_scale_weights(weights_dict, target_max=1000)

    session["weights_dict"] = weights_dict

    # Get classification from form (allow changing it during re-run)
    classification = request.form.get("classification")
    if classification:
        session["classification"] = classification

    # Get background from form (allow changing it during re-run)
    background_from_form = request.form.get("background")
    if background_from_form:
        # Same guard as in submission(), for the re-run form on the graph page.
        missing_leaves_json = _missing_background_file(background_from_form)
        if missing_leaves_json:
            return _render_background_unavailable(
                background_from_form,
                missing_leaves_json,
            )
        session["background"] = background_from_form

    # Get expand_background from form (allow changing it during re-run)
    if request.method == "POST" and "expand_background" in request.form:
        session["expand_background"] = bool(request.form.get("expand_background"))

    # Keep the correction selection stable across re-runs.
    previous_correction = session.get(
        "correction_method",
        {
            "bonferroni_correct": False,
            "benjamini_hochberg_correct": False,
        },
    )
    method = request.form.get("p_value_correction_method")
    if method:
        bonferroni_correct, benjamini_hochberg_correct = map_p_value_correction_method(
            method,
        )
    else:
        bonferroni_correct = previous_correction.get("bonferroni_correct", False)
        benjamini_hochberg_correct = previous_correction.get(
            "benjamini_hochberg_correct",
            False,
        )

    session["correction_method"] = {
        "bonferroni_correct": bonferroni_correct,
        "benjamini_hochberg_correct": benjamini_hochberg_correct,
    }

    background = session.get("background")
    looping_prune_method = request.form.get("looping_prune_method")

    if looping_prune_method == "no_loop_prune":
        # User chose no looping pruning, proceed to other pruning options
        pass
    elif looping_prune_method == "plain_enrich":
        if background == "full":
            # Use weighted analysis if weights are present, otherwise use standard analysis
            if weights_dict:
                results, pruned_G = (
                    run_weighted_enrichment_analysis_plain_enrich_pruning_strategy(
                        weights_dict,
                        classification=session.get("classification") or "structural",
                    )
                )
            else:
                results, pruned_G = (
                    run_enrichment_analysis_plain_enrich_pruning_strategy(
                        studyset_list,
                        classification=session.get("classification") or "structural",
                    )
                )
            # Save JSON representation of pruned_G in session for graph visualization
            graph_json_file = f"website/static/data/graph_{session['session_id']}.json"
            graph_to_cytospace_json(pruned_G, graph_json_file, results)
            session["graph_file"] = f"graph_{session['session_id']}.json"

            session["pruning"] = {
                "method": "plain_enrich",
            }

            session["unresolved_smiles"] = unresolved_smiles
            session["ambiguous_smiles_matches"] = ambiguous_smiles_matches

            return render_template(
                "results.html",
                results=results,
                graph_json_file=graph_json_file,
                unresolved_smiles=unresolved_smiles,
                ambiguous_smiles_matches=ambiguous_smiles_matches,
                smiles_option=session.get("smiles_option"),
                expand_background=session.get("expand_background", True),
                background=session.get("background", "full"),
            )

        elif background in NARROW_BACKGROUND_LEAVES:
            if weights_dict:
                (
                    results,
                    pruned_G,
                    leaves_to_expand_background,
                    parents_to_expand_background,
                ) = run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy(
                    weights_dict,
                    classification=session.get("classification") or "structural",
                    narrow_background_leaves_json=background,
                    expand_background=session.get("expand_background", True),
                )
            else:
                (
                    results,
                    pruned_G,
                    leaves_to_expand_background,
                    parents_to_expand_background,
                ) = run_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy(
                    studyset_list,
                    classification=session.get("classification") or "structural",
                    narrow_background_leaves_json=background,
                    expand_background=session.get("expand_background", True),
                )

            graph_json_file = f"website/static/data/graph_{session['session_id']}.json"
            graph_to_cytospace_json(pruned_G, graph_json_file, results)
            session["graph_file"] = f"graph_{session['session_id']}.json"

            session["pruning"] = {
                "method": "plain_enrich",
            }

            session["unresolved_smiles"] = unresolved_smiles
            session["ambiguous_smiles_matches"] = ambiguous_smiles_matches

            return render_template(
                "results.html",
                results=results,
                graph_json_file=graph_json_file,
                unresolved_smiles=unresolved_smiles,
                ambiguous_smiles_matches=ambiguous_smiles_matches,
                smiles_option=session.get("smiles_option"),
                leaves_to_expand_background=leaves_to_expand_background,
                parents_to_expand_background=parents_to_expand_background,
                expand_background=session.get("expand_background", True),
                background=session.get("background", "full"),
            )

    # PRUNING OPTIONS
    root_children_prune = request.form.get("root_children_prune") == "true"
    levels = int(request.form.get("levels", 2))

    linear_branch_prune = request.form.get("linear_branch_prune") == "true"
    linear_branch_n = int(request.form.get("linear_branch_n", 2))

    high_p_value_prune = request.form.get("high_p_value_prune") == "true"
    p_value_threshold = float(request.form.get("p_value_threshold", 0.05))

    zero_degree_prune = request.form.get("zero_degree_prune") == "true"

    background = session.get("background")

    # Initialize expanded leaves/parents tracking (for narrow backgrounds)
    leaves_to_expand_background = set()
    parents_to_expand_background = set()

    # Use weighted analysis if weights are present, otherwise use standard analysis
    if weights_dict:
        if background in NARROW_BACKGROUND_LEAVES:
            (
                results,
                pruned_G,
                leaves_to_expand_background,
                parents_to_expand_background,
            ) = run_weighted_narrow_background_enrichment_analysis(
                weights_dict,
                levels=levels,
                n=linear_branch_n,
                p_value_threshold=p_value_threshold,
                classification=session.get("classification") or "structural",
                root_children_prune=root_children_prune,
                linear_branch_prune=linear_branch_prune,
                high_p_value_prune=high_p_value_prune,
                zero_degree_prune=zero_degree_prune,
                bonferroni_correct=bonferroni_correct,
                benjamini_hochberg_correct=benjamini_hochberg_correct,
                narrow_background_leaves_json=background,
                expand_background=session.get("expand_background", True),
            )
        else:
            results, pruned_G = run_weighted_enrichment_analysis(
                weights_dict,
                levels=levels,
                n=linear_branch_n,
                p_value_threshold=p_value_threshold,
                classification=session.get("classification") or "structural",
                root_children_prune=root_children_prune,
                linear_branch_prune=linear_branch_prune,
                high_p_value_prune=high_p_value_prune,
                zero_degree_prune=zero_degree_prune,
                bonferroni_correct=bonferroni_correct,
                benjamini_hochberg_correct=benjamini_hochberg_correct,
            )
    else:
        if background in NARROW_BACKGROUND_LEAVES:
            (
                results,
                pruned_G,
                leaves_to_expand_background,
                parents_to_expand_background,
            ) = run_narrow_background_enrichment_analysis(
                studyset_list,
                bonferroni_correct=bonferroni_correct,
                benjamini_hochberg_correct=benjamini_hochberg_correct,
                root_children_prune=root_children_prune,
                levels=levels,
                linear_branch_prune=linear_branch_prune,
                n=linear_branch_n,
                high_p_value_prune=high_p_value_prune,
                p_value_threshold=p_value_threshold,
                zero_degree_prune=zero_degree_prune,
                classification=session.get("classification") or "structural",
                narrow_background_leaves_json=background,
                expand_background=session.get("expand_background", True),
            )
        else:
            results, pruned_G = run_enrichment_analysis(
                studyset_list,
                bonferroni_correct=bonferroni_correct,
                benjamini_hochberg_correct=benjamini_hochberg_correct,
                root_children_prune=root_children_prune,
                levels=levels,
                linear_branch_prune=linear_branch_prune,
                n=linear_branch_n,
                high_p_value_prune=high_p_value_prune,
                p_value_threshold=p_value_threshold,
                zero_degree_prune=zero_degree_prune,
                classification=session.get("classification") or "structural",
            )

    # Save JSON representation of pruned_G in session for graph visualization
    graph_json_file = f"website/static/data/graph_{session['session_id']}.json"
    graph_to_cytospace_json(pruned_G, graph_json_file, results)
    session["graph_file"] = f"graph_{session['session_id']}.json"

    session["pruning"] = {
        "method": "custom",
        "root_children_prune": root_children_prune,
        "levels": levels,
        "linear_branch_prune": linear_branch_prune,
        "linear_branch_n": linear_branch_n,
        "high_p_value_prune": high_p_value_prune,
        "p_value_threshold": p_value_threshold,
        "zero_degree_prune": zero_degree_prune,
    }

    session["unresolved_smiles"] = unresolved_smiles
    session["ambiguous_smiles_matches"] = ambiguous_smiles_matches

    return render_template(
        "results.html",
        results=results,
        graph_json_file=graph_json_file,
        unresolved_smiles=unresolved_smiles,
        ambiguous_smiles_matches=ambiguous_smiles_matches,
        smiles_option=session.get("smiles_option"),
        leaves_to_expand_background=leaves_to_expand_background,
        parents_to_expand_background=parents_to_expand_background,
        expand_background=session.get("expand_background", True),
        background=session.get("background", "full"),
    )


@app.route("/graph")
def graph():
    session_id = session.get("session_id")
    graph_file = f"graph_{session_id}.json" if session_id else "graph.json"
    pruning = session.get("pruning", {})
    correction_method = session.get("correction_method", {})
    classification = session.get("classification", "structural")
    background = session.get("background", "full")
    return render_template(
        "graph.html",
        graph_file=graph_file,
        pruning=pruning,
        correction_method=correction_method,
        classification=classification,
        background=background,
        smiles_option=session.get("smiles_option"),
        expand_background=session.get("expand_background", True),
    )


if __name__ == "__main__":
    app.run(debug=True)
