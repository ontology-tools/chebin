"""Export an enrichment graph as a self-contained interactive HTML page.

The website renders the graph through Flask, which package users do not have. This module
writes the same Cytoscape view to a single HTML file that opens in any browser, with no
server and no network access: Cytoscape is inlined from the vendored copy and the graph
data is embedded in the page.

The page is regenerated from the website template by
``tools/build_standalone_graph_template.py`` -- see that script for what differs.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

from chebin.calculations.visualitations_and_pruning import graph_to_cytoscape_dict

_PACKAGE_DIR = Path(__file__).resolve().parent
_TEMPLATE = _PACKAGE_DIR / "templates" / "graph_standalone.html"
_CYTOSCAPE_JS = _PACKAGE_DIR / "static" / "cytoscape.min.js"

_DATA_TOKEN = "/*__CHEBIN_GRAPH_DATA__*/"
_CYTOSCAPE_TOKEN = "/*__CHEBIN_CYTOSCAPE_JS__*/"


def _read_asset(path: Path, what: str) -> str:
    try:
        return path.read_text()
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{what} is missing from the chebin package at {path}. The installed package "
            f"looks incomplete -- reinstall it, or check that package data is included in "
            f"the build.",
        ) from None


def _embed_json(data: dict) -> str:
    """Serialise for embedding inside a <script> tag.

    ``</script>`` appearing inside a string literal would close the tag early, so the
    forward slash is escaped. The result is still valid JSON to `JSON.parse` and a valid
    JS literal.
    """
    return json.dumps(data).replace("</", "<\\/")


def export_graph_html(
    G,
    enrichment_results,
    output_file,
    include_untested_leaves=False,
):
    """Write the graph as a standalone interactive HTML page.

    The arguments mirror
    :func:`chebin.calculations.visualitations_and_pruning.graph_to_cytospace_json`.

    Args:
        G: The pruned graph returned by any of the ``run_*_enrichment_analysis``
            functions.
        enrichment_results: The results dict returned alongside `G`. Pass ``None`` if
            you don't have one -- the graph still renders, but no node carries a p-value.
        output_file: Path of the ``.html`` file to write.
        include_untested_leaves: Keep study-set leaf nodes, which are never tested and so
            never coloured. Excluded by default -- they are the large majority of nodes.

    Returns:
        str: The path written.

    Example:
        >>> results, graph = run_enrichment_analysis(["CHEBI:15377"])
        >>> export_graph_html(graph, results, "enrichment.html")
        'enrichment.html'
    """
    if enrichment_results is None:
        warnings.warn(
            "export_graph_html called without enrichment_results -- every node in "
            "the exported graph will have no p-value.",
            stacklevel=2,
        )
    data = graph_to_cytoscape_dict(
        G,
        enrichment_results=enrichment_results,
        include_untested_leaves=include_untested_leaves,
    )

    page = _read_asset(_TEMPLATE, "The standalone graph template")
    page = page.replace(
        _CYTOSCAPE_TOKEN,
        _read_asset(_CYTOSCAPE_JS, "The vendored Cytoscape bundle"),
    )
    page = page.replace(_DATA_TOKEN, _embed_json(data))

    output_path = Path(output_file)
    if output_path.parent != Path(""):
        output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(page)

    n_nodes = sum(1 for e in data["elements"] if "source" not in e.get("data", {}))
    print(f"Wrote standalone graph ({n_nodes} nodes) to {output_path}")
    return str(output_path)
