"""Every public enrichment entry point runs end to end against a real data folder.

The rest of the suite uses mock data and never calls these functions for real, so a
refactor can leave one of them referencing a name that no longer exists and every test
still passes -- which is exactly what happened when the duplicated data-file block was
replaced by `load_data_files()`: `run_narrow_background_enrichment_analysis` raised
`NameError: class_to_leaf_map_file` the first time anyone called it.

These tests need the generated data files, so they skip unless a data folder is
available. Point at one with CHEBIN_DATA_DIR, or run from a checkout with `data/` built.
"""

import os
from pathlib import Path

import pytest

from chebin.config import set_data_dir

REQUIRED_FILES = [
    "class_to_leaf_descendants_map.json",
    "class_to_all_roles_map.json",
    "roles_to_leaves_map.json",
    "removed_leaf_classes_with_smiles.csv",
    "removed_leaf_classes_to_ALL_parents_map.json",
    "chebi_parent_map.json",
    "chebi_id_to_name_map.json",
    "recon3d_leaves.json",
]

STUDY_SET = ["CHEBI:15377", "CHEBI:16236", "CHEBI:17234", "CHEBI:16947"]
WEIGHTS = {
    "CHEBI:15377": 2.5,
    "CHEBI:16236": 1.8,
    "CHEBI:17234": 3.1,
    "CHEBI:16947": 0.9,
}


def _find_data_dir():
    """A data folder with everything these tests need, or None."""
    candidates = []
    if os.environ.get("CHEBIN_DATA_DIR"):
        candidates.append(Path(os.environ["CHEBIN_DATA_DIR"]))
    candidates.append(Path.cwd() / "data")
    # tests/ -> binche2/data, so the suite works from any cwd in a checkout
    candidates.append(Path(__file__).resolve().parents[1] / "data")

    for candidate in candidates:
        if all((candidate / name).exists() for name in REQUIRED_FILES):
            return candidate
    return None


DATA_DIR = _find_data_dir()

pytestmark = pytest.mark.skipif(
    DATA_DIR is None,
    reason=(
        "needs a generated data folder; set CHEBIN_DATA_DIR or run create_all_files()"
    ),
)


@pytest.fixture(autouse=True)
def _point_at_data():
    set_data_dir(DATA_DIR)
    yield
    set_data_dir(None)


@pytest.fixture
def recon3d():
    """A narrow background that create_all_files always builds (unlike the HMDB one)."""
    return str(DATA_DIR / "recon3d_leaves.json")


def _assert_enrichment_shape(result, graph):
    assert set(result) >= {"study_set", "enrichment_results"}
    assert isinstance(result["enrichment_results"], dict)
    assert graph.number_of_nodes() > 0
    assert result["enrichment_results"], "no classes were tested"


def test_run_enrichment_analysis():
    from chebin.calculations.fishers_calculations import run_enrichment_analysis

    _assert_enrichment_shape(*run_enrichment_analysis(STUDY_SET))


def test_run_enrichment_analysis_plain_enrich():
    from chebin.calculations.fishers_calculations import (
        run_enrichment_analysis_plain_enrich_pruning_strategy as run,
    )

    _assert_enrichment_shape(*run(STUDY_SET))


def test_run_weighted_enrichment_analysis():
    from chebin.calculations.weighted_calculations import (
        run_weighted_enrichment_analysis,
    )

    _assert_enrichment_shape(*run_weighted_enrichment_analysis(WEIGHTS))


def test_run_weighted_narrow_background(recon3d):
    from chebin.calculations.weighted_calculations import (
        run_weighted_narrow_background_enrichment_analysis as run,
    )

    result, graph, *_ = run(WEIGHTS, narrow_background_leaves_json=recon3d)
    _assert_enrichment_shape(result, graph)


def test_run_weighted_narrow_background_plain_enrich(recon3d):
    from chebin.calculations.weighted_calculations import (
        run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy as run,
    )

    result, graph, *_ = run(WEIGHTS, narrow_background_leaves_json=recon3d)
    _assert_enrichment_shape(result, graph)


def test_run_narrow_background(recon3d):
    from chebin.preparing_data.wikidata.narrow_background_fishers import (
        run_narrow_background_enrichment_analysis as run,
    )

    result, graph, *_ = run(STUDY_SET, narrow_background_leaves_json=recon3d)
    _assert_enrichment_shape(result, graph)


def test_run_narrow_background_plain_enrich(recon3d):
    from chebin.preparing_data.wikidata.narrow_background_fishers import (
        run_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy as run,
    )

    result, graph, *_ = run(STUDY_SET, narrow_background_leaves_json=recon3d)
    _assert_enrichment_shape(result, graph)


def test_export_graph_html_is_self_contained(tmp_path):
    """The packaged graph export must not reach the network to render.

    Nothing may be *loaded* from the network -- no external script, stylesheet or
    image. The per-node ChEBI links are not that: they are inert until clicked, and
    each URL is assembled at runtime from the node's id, so no external href is
    written into the file either. See test_export_graph_html_links_nodes_to_chebi.
    """
    from chebin.calculations.fishers_calculations import run_enrichment_analysis
    from chebin.visualization import export_graph_html

    result, graph = run_enrichment_analysis(STUDY_SET)
    out = tmp_path / "graph.html"
    export_graph_html(graph, result, str(out))

    page = out.read_text()
    assert 'src="http' not in page
    assert 'href="http' not in page
    assert "url_for" not in page
    assert "/*__CHEBIN_GRAPH_DATA__*/" not in page
    assert "/*__CHEBIN_CYTOSCAPE_JS__*/" not in page
    assert "Cytoscape Consortium" in page  # the vendored bundle is inlined


def test_export_graph_html_links_nodes_to_chebi(tmp_path):
    """The exported page can take a node to its ChEBI entry.

    The link is built in JS from each node's id, so what the file can show is the
    machinery: the entry-page base URL and the two ways in (the tooltip title, and
    double-click).
    """
    from chebin.calculations.fishers_calculations import run_enrichment_analysis
    from chebin.visualization import export_graph_html

    result, graph = run_enrichment_analysis(STUDY_SET)
    out = tmp_path / "graph.html"
    export_graph_html(graph, result, str(out))

    page = out.read_text()
    assert "https://www.ebi.ac.uk/chebi/" in page
    assert "titleHtml" in page  # the tooltip title, with its id linked
    assert "dbltap" in page  # the double-click shortcut

    # The ids the link is built from survive into the inlined graph data.
    assert '"chebi_id": "CHEBI:' in page
