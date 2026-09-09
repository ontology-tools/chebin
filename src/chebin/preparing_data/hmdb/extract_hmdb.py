"""Extract HMDB metabolites into a flat tabular file.

The HMDB XML is large, so this script streams metabolite entries one at a time
instead of loading the whole file into memory.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Iterator
from pathlib import Path

HMDB_NAMESPACE = "http://www.hmdb.ca"
NS = {"hmdb": HMDB_NAMESPACE}
METABOLITE_TAG = f"{{{HMDB_NAMESPACE}}}metabolite"


def _clean_text(value: str | None) -> str:
    """Return stripped text or an empty string."""

    if value is None:
        return ""
    return value.strip()


def _iter_metabolites(
    xml_file: str,
    allow_truncated: bool = True,
) -> Iterator[ET.Element]:
    """Yield metabolite elements from the HMDB XML file."""

    context = ET.iterparse(xml_file, events=("end",))
    try:
        for _, element in context:
            if element.tag == METABOLITE_TAG:
                yield element
                element.clear()
    except ET.ParseError as error:
        if allow_truncated:
            line = "?"
            column = "?"
            if hasattr(error, "position"):
                line = str(error.position[0])
                column = str(error.position[1])
            print(
                f"Warning: XML parsing stopped near line {line}, column {column}. "
                "The input appears truncated; writing rows parsed so far.",
                file=sys.stderr,
            )
            return
        raise


def _get_text(element: ET.Element, path: str) -> str:
    """Get text from a namespaced child element."""

    return _clean_text(element.findtext(path, default="", namespaces=NS))


def _get_chebi_ids(metabolite: ET.Element) -> str:
    """Collect all non-empty ChEBI IDs for a metabolite."""

    chebi_ids = []
    for chebi_element in metabolite.findall(".//hmdb:chebi_id", NS):
        chebi_id = _clean_text(chebi_element.text)
        if chebi_id:
            chebi_ids.append(chebi_id)

    # Keep one row per metabolite, but preserve all available ChEBI IDs.
    return ";".join(dict.fromkeys(chebi_ids))


def _get_inchikey(metabolite: ET.Element) -> str:
    """Return the first non-empty InChIKey for a metabolite, or empty string."""
    # inchikey appears as a direct child <inchikey> under metabolite in HMDB XML
    return _clean_text(metabolite.findtext("hmdb:inchikey", default="", namespaces=NS))


def extract_hmdb_metabolites(
    xml_file: str,
    allow_truncated: bool = True,
) -> Iterator[dict[str, str]]:
    """Extract flat metabolite rows from the HMDB XML file."""

    for metabolite in _iter_metabolites(xml_file, allow_truncated=allow_truncated):
        yield {
            "hmdb_id": _get_text(metabolite, "hmdb:accession"),
            "name": _get_text(metabolite, "hmdb:name"),
            "chebi_id": _get_chebi_ids(metabolite),
            "smiles": _get_text(metabolite, "hmdb:smiles"),
            "inchikey": _get_inchikey(metabolite),
            "status": _get_text(metabolite, "hmdb:status"),
        }


def write_hmdb_table(
    rows: Iterable[dict[str, str]],
    output_file: str,
    format: str = "tsv",
) -> int:
    """Write extracted rows to a TSV or CSV file and return the row count."""

    output_path = Path(output_file)
    if format not in {"tsv", "csv"}:
        raise ValueError("format must be 'tsv' or 'csv'")

    delimiter = "\t" if format == "tsv" else ","
    output_path.parent.mkdir(parents=True, exist_ok=True)

    row_count = 0

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["hmdb_id", "name", "chebi_id", "smiles", "inchikey", "status"],
            delimiter=delimiter,
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            row_count += 1

    return row_count


def extract_hmdb_to_file(
    xml_file: str,
    output_file: str,
    format: str = "tsv",
    allow_truncated: bool = True,
) -> None:
    """Extract HMDB metabolites and write them directly to disk."""

    start_time = time.perf_counter()
    print(f"Starting HMDB extraction from {xml_file}...")
    row_count = write_hmdb_table(
        extract_hmdb_metabolites(xml_file, allow_truncated=allow_truncated),
        output_file,
        format=format,
    )
    elapsed_seconds = time.perf_counter() - start_time
    print(
        f"Finished writing {row_count} rows to {output_file} in "
        f"{elapsed_seconds:.1f} seconds ({elapsed_seconds / 60:.1f} minutes).",
    )


def _format_parse_error(xml_file: str, error: ET.ParseError) -> str:
    """Create a concise, actionable XML parse error message."""

    line = "?"
    column = "?"
    if hasattr(error, "position"):
        line = str(error.position[0])
        column = str(error.position[1])

    return (
        f"Failed to parse XML in {xml_file} near line {line}, column {column}. "
        "The HMDB file is likely truncated or malformed. "
        "Please re-download or validate the XML, then rerun the script."
    )


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract HMDB metabolites into a flat TSV or CSV file.",
    )
    parser.add_argument(
        "xml_file",
        nargs="?",
        default="data/hmdb_metabolites.xml",
        help="Path to the HMDB metabolites XML file.",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        default="data/hmdb_metabolites_extract.tsv",
        help="Output file path.",
    )
    parser.add_argument(
        "--format",
        choices=("tsv", "csv"),
        default="tsv",
        help="Output format. TSV is the default because it is safer for text fields.",
    )
    parser.add_argument(
        "--strict-xml",
        action="store_true",
        help="Fail on malformed/truncated XML instead of writing partial output.",
    )
    return parser


def create_hmdb_output_file() -> None:
    args = build_argument_parser().parse_args()
    try:
        extract_hmdb_to_file(
            args.xml_file,
            args.output_file,
            format=args.format,
            allow_truncated=not args.strict_xml,
        )
    except ET.ParseError as error:
        print(_format_parse_error(args.xml_file, error), file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    create_hmdb_output_file()
