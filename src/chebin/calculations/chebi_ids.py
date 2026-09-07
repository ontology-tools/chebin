"""Recognising and canonicalising the ChEBI identifiers people actually type.

Users paste ChEBI IDs in whichever form they copied them from: ``CHEBI:17079``,
``chebi:17079``, ``ChEBI:17079``, ``CHEBI ID: 17079``, ``CHEBI_17079``, a full OBO
IRI, or just ``17079``. Everything downstream expects a single form -- the OBO IRI
``http://purl.obolibrary.org/obo/CHEBI_17079`` -- so all the variants are folded
into it here, in one place shared by the website's study-set parser,
:func:`chebin.calculations.fishers_calculations.normalize_id` (which every
analysis entry point runs its study set through) and the SMILES/ChEBI-ID
discrimination in :mod:`chebin.calculations.smiles_lookup`.

A bare number counts as a ChEBI ID: nothing else in a study set is a plain
integer, since a SMILES is never all digits.
"""

from __future__ import annotations

import re

#: Namespace every ChEBI class IRI lives in.
CHEBI_IRI_PREFIX = "http://purl.obolibrary.org/obo/"

#: The prefix, in the ways people write it: ``CHEBI``/``ChEBI``/``chebi``, an
#: optional ``ID``, and any of the usual separators (or none at all, as in
#: ``CHEBI17079``). Shared by the whole-token and in-line patterns below.
_PREFIX_PATTERN = r"chebi\s*(?:id)?\s*[:_\s-]?\s*"

#: A whole token that is a ChEBI ID -- the prefix above, or no prefix at all.
_CHEBI_TOKEN_RE = re.compile(rf"^(?:{_PREFIX_PATTERN})?(\d+)$", re.IGNORECASE)

#: The same prefixes matched anywhere in a line, so that the spellings containing
#: whitespace (``CHEBI ID: 17079``) can be folded together *before* the line is
#: split into tokens. The lookbehind keeps IRIs (``.../CHEBI_17079``) out of it:
#: they are already canonical, and parse_chebi_number() handles them directly.
_CHEBI_INLINE_RE = re.compile(
    rf"(?<![/\w]){_PREFIX_PATTERN}(\d+)(?!\w)",
    re.IGNORECASE,
)


def parse_chebi_number(value: object) -> str | None:
    """The numeric part of a ChEBI identifier, whichever way it was written.

    Args:
        value: A candidate identifier, e.g. ``CHEBI:17079``, ``chebi 17079``,
            ``CHEBI ID: 17079``, ``CHEBI_17079``, ``17079``, or an IRI such as
            ``http://purl.obolibrary.org/obo/CHEBI_17079``.

    Returns:
        str | None: The ChEBI number without leading zeros (``"17079"``), or None
            if the value isn't a ChEBI identifier in any recognised form.
    """
    if value is None:
        return None

    text = str(value).strip().strip("\"'").strip()
    if not text:
        return None

    # For an IRI, only the last path segment can be the identifier.
    if "://" in text:
        text = text.rsplit("/", 1)[-1]

    match = _CHEBI_TOKEN_RE.match(text)
    if match is None:
        return None
    return str(int(match.group(1)))


def looks_like_chebi_id(value: object) -> bool:
    """Whether ``value`` is a ChEBI identifier rather than something else (a SMILES)."""
    return parse_chebi_number(value) is not None


def to_chebi_curie(value: object) -> str | None:
    """``value`` as ``CHEBI:17079``, or None if it isn't a ChEBI identifier."""
    number = parse_chebi_number(value)
    return None if number is None else f"CHEBI:{number}"


def to_chebi_iri(value: object) -> str | None:
    """``value`` as a full OBO IRI, or None if it isn't a ChEBI identifier."""
    number = parse_chebi_number(value)
    return None if number is None else f"{CHEBI_IRI_PREFIX}CHEBI_{number}"


def collapse_chebi_prefixes(text: str) -> str:
    """Rewrite every loosely written ChEBI ID in ``text`` as ``CHEBI:17079``.

    Study-set lines are split on whitespace and commas, which would tear
    ``CHEBI ID: 17079`` into three tokens that are individually meaningless.
    Running a line through this first makes each ChEBI ID a single token, so the
    split (and the id/weight columns it produces) still line up.

    IRIs are left alone -- they are already unambiguous, and a ``CHEBI:`` in the
    middle of one would not survive being rewritten.
    """
    return _CHEBI_INLINE_RE.sub(lambda match: f"CHEBI:{int(match.group(1))}", text)
