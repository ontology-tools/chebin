"""Regenerate the package's standalone graph template from the website's graph page.

`chebin.visualization.export_graph_html` ships a self-contained copy of the Cytoscape
view so package users can look at a graph without running Flask. The website's
`graph.html` stays the source of truth for the visualisation logic; this script re-derives
the package template from it, applying the differences a serverless page needs:

* the graph JSON is inlined instead of fetched from a Flask `static` URL,
* Cytoscape is inlined from the vendored copy instead of loaded from a CDN,
* PDF export is dropped -- it was jsPDF's only use, and PNG/JPEG come from Cytoscape
  itself,
* the Bootstrap navbar dropdowns become a plain toolbar, and the Flask re-run form and
  its modal are omitted.

Run after changing the visualisation in `website/templates/graph.html`:

    python tools/build_standalone_graph_template.py

The generated file is committed; nothing regenerates it at runtime.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

BINCHE2 = Path(__file__).resolve().parents[1]
GRAPH_HTML = BINCHE2 / "website" / "templates" / "graph.html"
STYLES_CSS = BINCHE2 / "website" / "static" / "css" / "styles.css"
OUTPUT = BINCHE2 / "src" / "chebin" / "templates" / "graph_standalone.html"

DATA_TOKEN = "/*__CHEBIN_GRAPH_DATA__*/"
CYTOSCAPE_TOKEN = "/*__CHEBIN_CYTOSCAPE_JS__*/"

# --- the pieces lifted out of graph.html ------------------------------------------

BODY_START = '<div class="graph-wrapper">'
BODY_END = "</div>\n</div>"

# --- substitutions applied to the lifted JS ---------------------------------------

FETCH_OLD = """    fetch("{{url_for('static', filename='data/' ~ graph_file)}}")
    .then(response => response.json())
    .then(data => {"""
FETCH_NEW = """    Promise.resolve(window.__CHEBIN_GRAPH_DATA__)
    .then(data => {"""

# PDF export removal. jsPDF's only job was wrapping an already-generated PNG.
PDF_SUBSTITUTIONS = [
    (
        "            const isPdf = format === 'pdf';\n",
        "",
    ),
    (
        """            const savePdf = (dataUrl, width, height) => {
                const { jsPDF } = window.jspdf;
                const pdf = new jsPDF({
                    orientation: width >= height ? 'landscape' : 'portrait',
                    unit: 'px',
                    format: [width, height],
                });
                pdf.addImage(dataUrl, 'PNG', 0, 0, width, height);
                pdf.save('graph.pdf');
            };

""",
        "",
    ),
    (
        """            if (!includeLegend) {
                if (isPdf) {
                    const img = new Image();
                    img.onload = () => savePdf(rawData, img.width, img.height);
                    img.src = rawData;
                } else {
                    downloadDataUrl(rawData);
                }
                return;
            }""",
        """            if (!includeLegend) {
                downloadDataUrl(rawData);
                return;
            }""",
    ),
    (
        """                if (isPdf) {
                    savePdf(canvas.toDataURL('image/png'), canvas.width, canvas.height);
                    return;
                }

""",
        "",
    ),
    (
        "                if (format === 'png' || format === 'jpeg' || format === 'pdf') {",
        "                if (format === 'png' || format === 'jpeg') {",
    ),
    (
        """        // PDF embeds the same raster image as a single page sized to match it.
""",
        "",
    ),
    (
        """            // PDF always embeds a lossless PNG raster; png/jpeg use their own native format.
""",
        "",
    ),
]

# --- CSS lifted from the website stylesheet ---------------------------------------

CSS_SELECTORS = [
    ".graph-wrapper",
    "#cy",
    "#node-tooltip",
    "#node-tooltip a",
    "#color-legend",
    "#color-legend h4",
    ".legend-container",
    ".legend-gradient",
    ".legend-labels",
    "#p-filter",
    ".p-filter-handle",
    "#p-filter > *:not(.p-filter-handle)",
    "#p-filter label",
    '#p-filter input[type="range"]',
    ".p-filter-ends",
    "#p-filter .form-check-label",
    "#p-filter p",
    # The on-canvas controls live inside .graph-wrapper, so they are lifted with
    # the body and need their styling lifted too.
    "#graph-tools",
    ".tools-handle",
    "#graph-tools button",
    "#graph-tools button:hover:enabled",
    "#graph-tools button:disabled",
    '#graph-tools button[aria-pressed="true"]',
    "#graph-tools button:focus-visible",
    ".tools-status",
]

# Standalone-only chrome: replaces base.html's Bootstrap and the navbar dropdowns.
STANDALONE_CSS = """
/* --- standalone page chrome (replaces base.html + Bootstrap) --- */
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
  color: #212529;
  background: #fff;
}
.text-muted { color: #6c757d; }
.form-check { display: flex; align-items: center; gap: 6px; }
.form-check-input { margin: 0; }
.form-check-label { margin: 0; }

#toolbar {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 18px;
  padding: 10px 20px;
  border-bottom: 1px solid #dee2e6;
  background: #f8f9fa;
  font-size: 13px;
}
#toolbar .group { display: flex; align-items: center; gap: 6px; }
#toolbar .group > strong { color: #555; margin-right: 2px; }
#toolbar button {
  font: inherit;
  padding: 3px 9px;
  border: 1px solid #ccc;
  border-radius: 4px;
  background: #fff;
  cursor: pointer;
}
#toolbar button:hover { background: #eee; }
/* Toggle state, which Bootstrap supplies on the website. */
#toolbar button[aria-pressed="true"] { background: #0f808d; border-color: #0f808d; color: #fff; }
#toolbar button:focus-visible { outline: 2px solid #0f808d; outline-offset: 2px; }
"""

# The navbar dropdowns need Bootstrap JS; a plain toolbar carries the same hooks
# (.export-option / .layout-option / .show-hide-option and #exportIncludeLegend)
# that the lifted JS already binds to, so the JS itself needs no changes.
TOOLBAR_HTML = """<div id="toolbar">
  <div class="group">
    <strong>Export</strong>
    <label class="form-check">
      <input class="form-check-input" type="checkbox" id="exportIncludeLegend" checked>
      <span class="form-check-label">Include scale</span>
    </label>
    <button type="button" class="export-option" data-format="png">PNG</button>
    <button type="button" class="export-option" data-format="jpeg">JPEG</button>
  </div>
  <div class="group">
    <strong>Layout</strong>
    <button type="button" class="layout-option" data-layout="breadthfirst">Breadthfirst</button>
    <button type="button" class="layout-option" data-layout="cose">Cose</button>
    <button type="button" class="layout-option" data-layout="grid">Grid</button>
    <button type="button" class="layout-option" data-layout="circle">Circle</button>
  </div>
  <div class="group">
    <strong>Show/Hide</strong>
    <button type="button" class="show-hide-option" data-action="show-significant-with-paths">Most relevant paths</button>
    <button type="button" class="show-hide-option" data-action="show-all">Show all</button>
    <button type="button" class="show-hide-option" data-action="hide-unselected">Selected only</button>
    <button type="button" class="show-hide-option" data-action="hide-selected">Hide selected</button>
    <button type="button" class="show-hide-option" data-action="hide-non-significant">Hide non-significant</button>
    <button type="button" class="show-hide-option" data-action="hide-labels" aria-pressed="false"><span class="control-text">Hide labels</span></button>
  </div>
</div>"""
# Note: "Remove selected" and "Reset graph" are deliberately absent here. They live
# in the #graph-tools panel inside .graph-wrapper, which this page lifts verbatim --
# repeating them would give the standalone page two buttons per action.


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def extract_css_rules(css_text):
    """Pull the named top-level rules out of the website stylesheet, in order."""
    out = []
    for selector in CSS_SELECTORS:
        pattern = re.compile(
            r"^" + re.escape(selector) + r"\s*\{.*?^\}",
            re.MULTILINE | re.DOTALL,
        )
        match = pattern.search(css_text)
        if not match:
            fail(f"CSS rule not found in styles.css: {selector}")
        out.append(match.group(0))
    return "\n\n".join(out)


def extract_body(html):
    start = html.index(BODY_START)
    end = html.index(BODY_END, start) + len(BODY_END)
    return html[start:end]


def extract_js(html):
    """The main <script> block: the one containing the graph fetch."""
    blocks = re.findall(r"<script>(.*?)</script>", html, re.DOTALL)
    main = [b for b in blocks if "fetch(" in b and "cytoscape(" in b]
    if len(main) != 1:
        fail(f"expected exactly 1 main <script> block, found {len(main)}")
    return main[0]


def apply_substitutions(js):
    for old, new in [
        (FETCH_OLD, FETCH_NEW),
        *PDF_SUBSTITUTIONS,
    ]:
        if js.count(old) != 1:
            fail(
                f"expected exactly 1 occurrence of:\n{old[:160]}\n(found {js.count(old)})",
            )
        js = js.replace(old, new)
    return js


def main():
    if not GRAPH_HTML.exists():
        fail(f"{GRAPH_HTML} not found")

    html = GRAPH_HTML.read_text()
    css = extract_css_rules(STYLES_CSS.read_text())
    body = extract_body(html)
    js = apply_substitutions(extract_js(html))

    # Calls, not mentions: the tooltip code explains at length why it does not use
    # cytoscape-qtip, and that prose is not a dependency on it.
    for banned in ("{{", "{%", "url_for", ".qtip(", "jspdf.", "jsPDF"):
        if banned in js or banned in body:
            fail(
                f"'{banned}' still present in the lifted markup -- it would not work offline",
            )

    page = f"""<!DOCTYPE html>
<!--
  Generated by tools/build_standalone_graph_template.py from
  website/templates/graph.html. Do not edit by hand; edit the website template
  and regenerate, or the two will drift.

  Self-contained: no network access required.
-->
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ChEBI-N enrichment graph</title>
<style>
{css}
{STANDALONE_CSS}</style>
<script>{CYTOSCAPE_TOKEN}</script>
<script>window.__CHEBIN_GRAPH_DATA__ = {DATA_TOKEN};</script>
</head>
<body>

{TOOLBAR_HTML}

{body}

<script>{js}</script>

</body>
</html>
"""

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(page)
    print(f"wrote {OUTPUT.relative_to(BINCHE2)}  ({len(page):,} bytes)")
    print(f"  css rules lifted : {len(CSS_SELECTORS)}")
    print(f"  js lifted        : {len(js):,} bytes")


if __name__ == "__main__":
    main()
