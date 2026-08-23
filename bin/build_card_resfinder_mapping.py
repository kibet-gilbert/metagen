#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_mapping.py -- CARD <-> ResFinder acquired-resistance-gene name reconciliation.

WHAT THIS DOES
--------------
CARD (the Comprehensive Antibiotic Resistance Database; McArthur et al. 2013,
Alcock et al. 2020/2023) and ResFinder (Zankari et al. 2012; Bortolaia et al. 2020)
both catalogue antimicrobial-resistance genes, but under different naming
conventions and with different scope. This script builds a merged dictionary
matching CARD Antibiotic Resistance Ontology (ARO) terms -- as exported in a CARD
RGI "wildcard" prevalence report -- to ResFinder DB acquired-gene names, using a
series of literature-grounded, deterministic rules rather than free-text fuzzy
matching:

  1. exact / case-insensitive string match
  2. van-cluster structural rule:
       CARD "vanH gene in vanA cluster"  <->  ResFinder "VanH-A"
  3. punctuation/case normalization (apostrophe COUNT is preserved, because it is
     semantically load-bearing for aminoglycoside-modifying-enzyme regiospecificity,
     e.g. APH(3')-Ia vs APH(3'')-Ia are different enzymes -- Shaw et al. 1993)
  4. ResFinder 'bla' prefix stripping:
       CARD "CTX-M-58"  <->  ResFinder "blaCTX-M-58"
  5. alias-bridge matching via "Alternative name" annotations mined from
     ResFinder's phenotypes.txt (recovers historical aminoglycoside gene names,
     e.g. ResFinder "strA" = CARD "APH(3'')-Ib" -- Ramirez & Tolmasky 2010)

A data-driven drug-class consistency check (learned from already-unambiguous
matches, not hand-typed) then re-ranks candidates to catch same-spelling /
different-gene collisions across unrelated resistance mechanisms, before a small,
explicitly logged table of literature-grounded manual overrides is applied.

Only CARD's "protein homolog model" entries (horizontally acquired genes) are
expected to find a ResFinder counterpart; "protein variant model" (chromosomal
point mutations), "rRNA gene variant model", and "protein overexpression model"
entries are categorically out of ResFinder's acquired-gene scope (that territory
belongs to ResFinder's companion tool, PointFinder -- Zankari et al. 2017) and are
reported as such rather than forced into a match.

For full methodology, validation results, and reference list see documentation in: /docs/.


INPUT FILES (as distributed alongside ResFinder DB / exported from CARD RGI)
-----------------------------------------------------------------------------
  card_prevalence.txt      CARD prevalence export. Tab-separated columns:
                            ARO Accession, Name, Model ID, Model Type, Pathogen,
                            NCBI Plasmid, NCBI WGS, NCBI Chromosome,
                            NCBI Genomic Island, Criteria, ARO Categories.
  notes.txt                ResFinder DB canonical acquired-gene list. Lines of
                            the form "gene_name:description:" grouped under
                            "#<Class>:" header lines.
  phenotypes.txt            ResFinder DB per-reference-sequence phenotype table.
                            Tab-separated columns include Gene_accession no.,
                            Class, Phenotype, PMID, Mechanism of resistance,
                            Notes, Required_gene. The free-text Notes column is
                            mined here for "Alternative name X, Y" annotations.

OUTPUT FILES (written to --outputdir)
--------------------------------------
  <prefix>.csv              One row per CARD ARO term: match status, matched
                             ResFinder gene (if any), match tier/confidence,
                             documented synonyms, excluded class-mismatched
                             candidates, and CARD NCBI-prevalence context.
  <prefix>.json              Same information as a keyed dictionary (CARD ARO
                             Accession -> details), plus a reverse ResFinder ->
                             CARD index.
  ambiguous_cases_review.csv  Every CARD entry that had more than one raw
                             candidate before re-ranking, with full provenance,
                             for manual audit.

USAGE EXAMPLES
---------------
  # Defaults (reads from /mnt/user-data/uploads, writes to /mnt/user-data/outputs)
  python3 build_mapping.py

  # Explicit directories
  python3 build_mapping.py --inputdir ./data --outputdir ./results

  # Point at individually renamed/relocated files
  python3 build_mapping.py \\
      --card-prevalence-file /data/card_v3.9/prevalence.tsv \\
      --resfinder-notes-file /data/resfinder_db/notes.txt \\
      --resfinder-phenotypes-file /data/resfinder_db/phenotypes.txt \\
      --outputdir ./results

  # Inspect the effect of the safety-net rules by disabling them
  python3 build_mapping.py --no-category-check --no-manual-overrides -v

Author: generated in an interactive CARD/ResFinder reconciliation session.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

__version__ = "1.1.0"

# ============================================================================
# Constants
# ============================================================================

#: Default filenames looked for under --inputdir when a specific --*-file
#: argument is not supplied on the command line.
DEFAULT_FILENAMES = {
    "card_prevalence": "card_prevalence.txt",
    "resfinder_notes": "notes.txt",
    "resfinder_phenotypes": "phenotypes.txt",
}

#: Match tiers in descending order of confidence. Lower index = more trusted.
#: See the module docstring / METHODOLOGY.md for what each tier means.
TIER_PRIORITY = ["exact", "case_insensitive", "van_structural", "loose", "debla", "alias_bridge"]

#: Confidence label surfaced in the output files for each match tier.
CONFIDENCE_BY_TIER = {
    "exact": "High",
    "case_insensitive": "High",
    "loose": "High",
    "debla": "High",
    "van_structural": "High",
    "alias_bridge": "Medium",
}

#: Literature-grounded manual curation overrides, applied after the automated
#: tiers and the drug-class consistency check. Keyed by CARD ARO Accession ->
#: the ResFinder gene name that should be promoted to "primary" match.
#:
#: ARO:3002639 CARD 'APH(3'')-Ib' ties between two alias_bridge candidates,
#: 'strA' and 'aph(6)-Ia', both reachable only because ResFinder's own
#: phenotypes.txt lists "strA" as an alternative name for BOTH families (a
#: shared third-party alias, not a direct family<->family link). The strA/strB
#: streptomycin-phosphotransferase operon is consistently documented as
#: APH(3'')-Ib / APH(6)-Id in the aminoglycoside-modifying-enzyme nomenclature
#: literature (Shaw et al. 1993 Microbiol Rev 57:138-63; Ramirez & Tolmasky
#: 2010 Drug Resist Updat 13:151-71) -- i.e. position-3'' phosphorylation,
#: matching CARD's own APH(3'') regiospecificity label. aph(6)-Ia (position-6
#: phosphorylation) is a different regiospecificity; its shared "strA" tag more
#: likely reflects pre-systematic-nomenclature literature overlap than true
#: identity, so 'strA' -- not 'aph(6)-Ia' -- is promoted to primary here.
MANUAL_OVERRIDES = {
    "ARO:3002639": "strA",
}


# ============================================================================
# Loading & parsing
# ============================================================================

def _require_columns(df: pd.DataFrame, columns: list[str], source: str) -> None:
    """Raise a clear error if expected columns are missing from a loaded file.

    Args:
        df: The loaded DataFrame to check.
        columns: Column names that must be present.
        source: Human-readable file path/description, used in the error message.

    Raises:
        SystemExit: if any required column is missing.
    """
    missing = [c for c in columns if c not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: '{source}' is missing expected column(s): {missing}\n"
            f"       Found columns: {list(df.columns)}\n"
            f"       Is this really a CARD/ResFinder file, or has its format changed?"
        )


def load_card_prevalence(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load CARD's RGI wildcard prevalence export.

    Args:
        path: Path to card_prevalence.txt (tab-separated).

    Returns:
        A tuple ``(card, card_unique)`` where ``card`` is the full, possibly
        multi-row-per-ARO-term table (one row per pathogen x criteria
        combination) and ``card_unique`` is deduplicated to one row per
        ``ARO Accession``.
    """
    card = pd.read_csv(path, sep="\t")
    _require_columns(
        card,
        ["ARO Accession", "Name", "Model Type", "Pathogen", "NCBI WGS", "Criteria", "ARO Categories"],
        path,
    )
    card_unique = card.drop_duplicates(subset=["ARO Accession"]).copy().reset_index(drop=True)
    return card, card_unique





def parse_resfinder_notes(path: str) -> pd.DataFrame:
    """Parse ResFinder DB's notes.txt into a tidy gene-name table.

    The file format is a series of ``gene_name:description:`` lines grouped
    under ``#<Class>:`` header lines (e.g. ``#Aminoglycoside:``).

    Args:
        path: Path to notes.txt.

    Returns:
        DataFrame with columns ``gene_name``, ``notes_class``, ``description``,
        one row per ResFinder canonical acquired-gene name.
    """
    with open(path, encoding="utf-8") as fh:
        lines = [line.strip("\r\n") for line in fh]

    records, current_class = [], None
    for line in lines:
        if not line.strip():
            continue
        if line.strip().startswith("#"):
            current_class = line.strip().lstrip("#").strip().rstrip(":").strip()
            continue
        parts = line.split(":")
        gene = parts[0].strip()
        desc = parts[1].strip() if len(parts) > 1 else ""
        records.append({"gene_name": gene, "notes_class": current_class, "description": desc})

    df = pd.DataFrame(records)
    if df.empty:
        sys.exit(f"ERROR: parsed zero gene entries from '{path}' -- is this really ResFinder's notes.txt?")
    return df


def _extract_alt_names(notes_value: object) -> list[str]:
    """Pull out 'Alternative name X, Y' cross-references from a free-text cell.

    ResFinder's phenotypes.txt 'Notes' column is an unstructured, semicolon-
    delimited bag of curator remarks. This looks for any semicolon-delimited
    chunk containing the phrase 'alternative name' (case-insensitive) and
    returns everything after it, comma-split into individual names.

    Args:
        notes_value: The raw cell value (may be NaN/non-string).

    Returns:
        List of alternative-name strings found in this cell (possibly empty).
    """
    if not isinstance(notes_value, str):
        return []
    found = []
    for chunk in notes_value.split(";"):
        idx = chunk.lower().find("alternative name")
        if idx == -1:
            continue
        rest = chunk[idx + len("alternative name"):]
        for name in rest.split(","):
            name = name.strip().strip(".").strip()
            if name:
                found.append(name)
    return found


def load_resfinder_phenotypes(path: str) -> pd.DataFrame:
    """Load ResFinder DB's phenotypes.txt and derive per-family alias data.

    Args:
        path: Path to phenotypes.txt (tab-separated).

    Returns:
        DataFrame with the original columns plus ``family`` (the gene family
        name recovered by stripping the trailing ``_<variant#>_<accession>``
        suffix from ``Gene_accession no.``) and ``alt_names`` (list of
        alternative names mined from the ``Notes`` column via
        :func:`_extract_alt_names`).
    """
    pheno = pd.read_csv(path, sep="\t")
    _require_columns(pheno, ["Gene_accession no.", "Notes"], path)
    pheno["family"] = pheno["Gene_accession no."].apply(
        lambda s: re.sub(r"_\d+_[A-Za-z0-9_.\-]+$", "", str(s))
    )
    pheno["alt_names"] = pheno["Notes"].apply(_extract_alt_names)
    return pheno


# ============================================================================
# Alias graph (union-find over phenotypes.txt "Alternative name" annotations)
# ============================================================================

class UnionFind:
    """Minimal union-find/disjoint-set structure used to group ResFinder gene
    families with the historical alternative names documented for them in
    phenotypes.txt (e.g. linking 'aph(3'')-Ib' with 'strA', or 'ant(2'')-Ia'
    with 'aadB')."""

    def __init__(self) -> None:
        self._parent: dict[str, str] = {}

    def find(self, x: str) -> str:
        """Return the representative (root) element of x's group, with path
        compression."""
        self._parent.setdefault(x, x)
        while self._parent[x] != x:
            self._parent[x] = self._parent[self._parent[x]]
            x = self._parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        """Merge the groups containing ``a`` and ``b``."""
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self._parent[ra] = rb

    def groups(self) -> dict[str, set[str]]:
        """Return ``{root: {members...}}`` for every group currently known."""
        out: dict[str, set[str]] = defaultdict(set)
        for node in self._parent:
            out[self.find(node)].add(node)
        return out

    def __contains__(self, x: str) -> bool:
        return x in self._parent


def build_alias_graph(pheno: pd.DataFrame) -> UnionFind:
    """Build the alias union-find structure from a loaded phenotypes.txt table.

    Args:
        pheno: DataFrame as returned by :func:`load_resfinder_phenotypes`
            (must have ``family`` and ``alt_names`` columns).

    Returns:
        A populated :class:`UnionFind` linking each gene family to every
        alternative name documented for it (and transitively, to any other
        family that shares one of those alternative names).
    """
    uf = UnionFind()
    for _, row in pheno.iterrows():
        family = row["family"]
        uf.union(family, family)
        for alt in row["alt_names"]:
            uf.union(family, alt)
    return uf


# ============================================================================
# Name normalization rules
# ============================================================================

def loose_key(name: str) -> str:
    """Normalize a gene name for case/punctuation-insensitive comparison.

    Lowercases and strips spaces, parentheses, and hyphens, but deliberately
    *keeps* apostrophes: apostrophe count is semantically load-bearing in
    aminoglycoside-modifying-enzyme nomenclature (regiospecificity -- e.g.
    APH(3')-Ia and APH(3'')-Ia are different enzymes; Shaw et al. 1993), so a
    naive alnum-only strip would wrongly collapse distinct genes together.

    Args:
        name: Raw gene/ARO-term name.

    Returns:
        The normalized key.
    """
    return re.sub(r"[^a-z0-9']", "", name.lower())


def debla(name: str) -> str:
    """Strip a leading ResFinder-style 'bla' beta-lactamase prefix, if present.

    CARD names beta-lactamase alleles without the 'bla' gene-symbol prefix
    (e.g. "CTX-M-58"), while ResFinder retains it ("blaCTX-M-58") -- both are
    standard conventions (Bush & Jacoby 2010). This makes the two comparable.

    Args:
        name: A (possibly bla-prefixed) gene name.

    Returns:
        ``name`` with a leading 'bla' removed (case-insensitive), or
        unchanged if it doesn't start with 'bla'.
    """
    if name.lower().startswith("bla") and len(name) > 3:
        return name[3:]
    return name


_VAN_RE_CARD = re.compile(r"^van([A-Za-z0-9]+) gene in van([A-Za-z0-9]+) cluster$", re.IGNORECASE)


def van_card_key(name: str) -> str | None:
    """Convert a CARD van-cluster descriptive name to a comparable key.

    CARD spells out auxiliary vanRSHAXYZ cluster genes descriptively, e.g.
    ``"vanH gene in vanA cluster"``; ResFinder uses the compact form
    ``"VanH-A"``. This produces the same normalized key for both, e.g. both
    become ``"vanha"``.

    Args:
        name: A CARD ARO term name.

    Returns:
        The normalized structural key if ``name`` matches the
        "vanX gene in vanY cluster" pattern, else ``None``.
    """
    m = _VAN_RE_CARD.match(name.strip())
    if not m:
        return None
    gene, cluster = m.group(1), m.group(2)
    return "van" + gene.lower() + cluster.lower()


# ============================================================================
# Matching
# ============================================================================

def build_resfinder_alias_table(notes_df: pd.DataFrame, alias_uf: UnionFind) -> pd.DataFrame:
    """Attach each ResFinder canonical gene name's full alias set.

    Args:
        notes_df: DataFrame from :func:`parse_resfinder_notes`.
        alias_uf: Alias graph from :func:`build_alias_graph`.

    Returns:
        DataFrame with columns ``canonical``, ``notes_class``, ``description``,
        ``aliases`` (sorted list including the canonical name itself plus any
        names linked to it via the alias graph).
    """
    groups = alias_uf.groups()
    records = []
    for _, row in notes_df.iterrows():
        canon = row["gene_name"]
        alias_set = {canon}
        if canon in alias_uf:
            alias_set |= groups[alias_uf.find(canon)]
        records.append(
            {
                "canonical": canon,
                "notes_class": row["notes_class"],
                "description": row["description"],
                "aliases": sorted(alias_set),
            }
        )
    return pd.DataFrame(records)


def build_resfinder_lookup_index(resfinder_df: pd.DataFrame) -> dict[str, list[tuple[str, str, str]]]:
    """Build a normalized-key -> candidate-match reverse index.

    Every alias of every ResFinder canonical gene name is indexed under both
    its direct normalized key and its bla-stripped normalized key, so a single
    dict lookup during matching covers the 'loose', 'debla', and
    'alias_bridge' tiers at once.

    Args:
        resfinder_df: DataFrame from :func:`build_resfinder_alias_table`.

    Returns:
        Mapping of ``normalized_key -> [(canonical_name, variant_tag, alias), ...]``
        where ``variant_tag`` is ``'direct'`` or ``'debla'``.
    """
    index: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for _, row in resfinder_df.iterrows():
        canon = row["canonical"]
        for alias in row["aliases"]:
            for variant, tag in [(alias, "direct"), (debla(alias), "debla")]:
                key = loose_key(variant)
                if key:
                    index[key].append((canon, tag, alias))
    return index


def match_card_to_resfinder(
    card_unique: pd.DataFrame,
    resfinder_df: pd.DataFrame,
    rf_index: dict[str, list[tuple[str, str, str]]],
) -> pd.DataFrame:
    """Run the full tiered matching pass for every CARD ARO term.

    Args:
        card_unique: One row per CARD ARO Accession (see :func:`load_card_prevalence`).
        resfinder_df: From :func:`build_resfinder_alias_table`.
        rf_index: From :func:`build_resfinder_lookup_index`.

    Returns:
        DataFrame with columns ``ARO_Accession``, ``CARD_Name``, ``Model_Type``,
        ``ARO_Categories``, ``n_candidates``, ``matches`` (list of dicts, each
        ``{'resfinder_name', 'tier_tag', 'notes_class'}``, sorted best-tier-first).
    """
    notes_class_map = dict(zip(resfinder_df["canonical"], resfinder_df["notes_class"]))
    canonical_set = set(resfinder_df["canonical"].values)
    canonical_lower_map: dict[str, list[str]] = defaultdict(list)
    for canon in canonical_set:
        canonical_lower_map[canon.lower()].append(canon)

    results = []
    for _, row in card_unique.iterrows():
        aro, name, mtype, cats = row["ARO Accession"], row["Name"], row["Model Type"], row["ARO Categories"]

        best_tier: dict[str, str] = {}

        def consider(canon: str, tier: str) -> None:
            if canon not in best_tier or TIER_PRIORITY.index(tier) < TIER_PRIORITY.index(best_tier[canon]):
                best_tier[canon] = tier

        # Tiers 1-2: exact / case-insensitive match against the canonical name itself.
        if name in canonical_set:
            consider(name, "exact")
        else:
            for canon in canonical_lower_map.get(name.lower(), []):
                consider(canon, "case_insensitive")

        # Tiers 3-4 & alias_bridge: normalized-key lookups (direct + debla variants).
        lk = loose_key(name)
        for canon, subtag, alias in rf_index.get(lk, []):
            if subtag == "direct":
                consider(canon, "loose" if alias == canon else "alias_bridge")
            elif subtag == "debla":
                consider(canon, "debla" if alias == canon else "alias_bridge")

        # van-cluster structural rule.
        vk = van_card_key(name)
        if vk:
            for canon, _subtag, _alias in rf_index.get(vk, []):
                consider(canon, "van_structural")

        resolved = [
            {"resfinder_name": canon, "tier_tag": tier, "notes_class": notes_class_map.get(canon)}
            for canon, tier in best_tier.items()
        ]
        resolved.sort(key=lambda m: TIER_PRIORITY.index(m["tier_tag"]))

        results.append(
            {
                "ARO_Accession": aro,
                "CARD_Name": name,
                "Model_Type": mtype,
                "ARO_Categories": cats,
                "n_candidates": len(resolved),
                "matches": resolved,
            }
        )
    return pd.DataFrame(results)


def apply_category_consistency_check(match_df: pd.DataFrame, enabled: bool = True) -> pd.DataFrame:
    """Re-rank each CARD entry's candidates using a data-driven drug-class check.

    Learns which CARD ``ARO Categories`` terms co-occur with which ResFinder
    ``notes_class`` values from already-unambiguous, high-confidence matches
    (exact / case_insensitive / debla / van_structural tiers with exactly one
    candidate), then demotes/excludes candidates that contradict that learned
    crosswalk. This is what catches same-spelling / different-gene collisions
    such as CARD's beta-lactamase 'FAR-1' string-matching ResFinder's unrelated
    fusidic-acid gene 'far1' (the correct match, 'blaFAR-1', is promoted
    instead), and correctly demotes CARD's 'SatA' to "no reliable match" rather
    than accepting its only, class-inconsistent alias-bridge candidate.

    Args:
        match_df: Output of :func:`match_card_to_resfinder`.
        enabled: If False, skip re-ranking entirely -- ``primary_match`` is
            just the best-tier raw candidate, with no class-consistency
            filtering. Useful for auditing what the safety net changes.

    Returns:
        ``match_df`` with three new columns added: ``primary_match`` (dict or
        None), ``alternate_matches`` (list of dicts), ``category_flagged``
        (list of dicts excluded for class inconsistency; always empty when
        ``enabled=False``).
    """
    if not enabled:
        primaries = [r["matches"][0] if r["matches"] else None for r in match_df.to_dict("records")]
        alternates = [r["matches"][1:] if r["matches"] else [] for r in match_df.to_dict("records")]
        match_df = match_df.copy()
        match_df["primary_match"] = primaries
        match_df["alternate_matches"] = alternates
        match_df["category_flagged"] = [[] for _ in range(len(match_df))]
        return match_df

    class_pair_counts: Counter = Counter()
    for _, r in match_df.iterrows():
        if r["n_candidates"] == 1 and r["matches"][0]["tier_tag"] in (
            "exact", "case_insensitive", "debla", "van_structural",
        ):
            m = r["matches"][0]
            for c in (c.strip() for c in r["ARO_Categories"].split(";")):
                class_pair_counts[(c, m["notes_class"])] += 1

    compatible: dict[str, set[str]] = defaultdict(set)
    for cat, nclass in class_pair_counts:
        compatible[nclass].add(cat)

    def category_consistent(aro_categories_str: str, notes_class: str) -> bool:
        cats = {c.strip() for c in aro_categories_str.split(";")}
        return bool(cats & compatible.get(notes_class, set()))

    def pick(row) -> tuple[dict | None, list[dict], list[dict]]:
        matches = row["matches"]
        if not matches:
            return None, [], []
        flagged = [m for m in matches if not category_consistent(row["ARO_Categories"], m["notes_class"])]
        consistent = [m for m in matches if m not in flagged]
        if not consistent:
            # Every candidate contradicts the learned crosswalk -- don't assert a
            # confident primary match; demote the whole entry to "no reliable match".
            return None, [], matches
        primary = consistent[0]
        alternates = [m for m in matches if m is not primary]
        return primary, alternates, flagged

    primaries, alternates_list, flagged_list = [], [], []
    for _, row in match_df.iterrows():
        p, a, f = pick(row)
        primaries.append(p)
        alternates_list.append(a)
        flagged_list.append(f)

    match_df = match_df.copy()
    match_df["primary_match"] = primaries
    match_df["alternate_matches"] = alternates_list
    match_df["category_flagged"] = flagged_list
    return match_df


def apply_manual_overrides(
    match_df: pd.DataFrame, overrides: dict[str, str], log, enabled: bool = True
) -> pd.DataFrame:
    """Apply the small, explicitly documented literature-grounded override table.

    Args:
        match_df: Output of :func:`apply_category_consistency_check`.
        overrides: ``{ARO_Accession: preferred_resfinder_name}``. See
            :data:`MANUAL_OVERRIDES` for the current table and its rationale.
        enabled: If False, this is a no-op (useful for auditing).
        log: Callable taking a single string, used to report each override
            applied (e.g. a logging function honoring --verbose/--quiet).

    Returns:
        ``match_df`` with ``primary_match``/``alternate_matches`` updated in
        place for any accession present in ``overrides``.
    """
    if not enabled:
        return match_df
    match_df = match_df.copy()
    applied = 0
    for idx, row in match_df.iterrows():
        target = overrides.get(row["ARO_Accession"])
        if target is None:
            continue
        candidate = next((m for m in row["matches"] if m["resfinder_name"] == target), None)
        if candidate is None:
            continue
        old_primary = row["primary_match"]
        match_df.at[idx, "primary_match"] = candidate
        match_df.at[idx, "alternate_matches"] = [m for m in row["matches"] if m is not candidate]
        applied += 1
        old_name = old_primary["resfinder_name"] if old_primary else None
        log(f"  manual override: {row['ARO_Accession']} '{row['CARD_Name']}': {old_name} -> {target}")
    log(f"Manual curation overrides applied: {applied}")
    return match_df


# ============================================================================
# Prevalence context & output assembly
# ============================================================================

def summarize_prevalence(card_full: pd.DataFrame) -> pd.DataFrame:
    """Collapse the full (multi-row-per-ARO-term) CARD table to one summary row
    per ARO Accession, for downstream comparative-prevalence work.

    Args:
        card_full: The un-deduplicated table from :func:`load_card_prevalence`
            (one row per pathogen x detection-criteria combination).

    Returns:
        DataFrame with columns ``ARO Accession``, ``n_pathogens_detected_in``,
        ``max_ncbi_wgs_pct``, ``top_pathogen_by_wgs_pct``, ``has_perfect_strict_hit``.
    """
    def _summarize(g: pd.DataFrame) -> pd.Series:
        idxmax = g["NCBI WGS"].idxmax()
        return pd.Series(
            {
                "n_pathogens_detected_in": g["Pathogen"].nunique(),
                "max_ncbi_wgs_pct": g["NCBI WGS"].max(),
                "top_pathogen_by_wgs_pct": g.loc[idxmax, "Pathogen"],
                "has_perfect_strict_hit": (g["Criteria"] == "perfect_strict").any(),
            }
        )

    return card_full.groupby("ARO Accession").apply(_summarize, include_groups=False).reset_index()


def _match_status(row: pd.Series) -> str:
    """Classify a matched/unmatched CARD entry for reporting.

    Distinguishes genuine "no reliable match found" (model type that, in
    principle, could have a ResFinder counterpart) from entries that are
    categorically out of ResFinder's acquired-gene scope by CARD model type.
    """
    if row["primary_match"] is not None:
        return "matched"
    if row["Model_Type"] == "protein variant model":
        return "out_of_scope_chromosomal_point_mutation"
    if row["Model_Type"] == "rRNA gene variant model":
        return "out_of_scope_rRNA_mutation"
    if row["Model_Type"] == "protein overexpression model":
        return "out_of_scope_regulatory_overexpression"
    return "no_reliable_match"


def build_merged_table(match_df: pd.DataFrame, prev_summary: pd.DataFrame) -> pd.DataFrame:
    """Assemble the final, flat, one-row-per-CARD-ARO-term output table.

    Args:
        match_df: Output of :func:`apply_manual_overrides`.
        prev_summary: Output of :func:`summarize_prevalence`.

    Returns:
        The merged DataFrame, sorted by (Match_Status, ARO_Accession), as
        written to the CSV/JSON deliverables.
    """
    rows = []
    for _, r in match_df.iterrows():
        prim = r["primary_match"]
        alts = r["alternate_matches"] or []
        flagged = r["category_flagged"] or []
        rows.append(
            {
                "ARO_Accession": r["ARO_Accession"],
                "CARD_Name": r["CARD_Name"],
                "Model_Type": r["Model_Type"],
                "ARO_Categories": r["ARO_Categories"],
                "Match_Status": _match_status(r),
                "ResFinder_Gene_Name": prim["resfinder_name"] if prim else "",
                "ResFinder_Class": prim["notes_class"] if prim else "",
                "Match_Tier": prim["tier_tag"] if prim else "",
                "Match_Confidence": CONFIDENCE_BY_TIER.get(prim["tier_tag"], "") if prim else "",
                "Alternate_ResFinder_Synonyms": "; ".join(a["resfinder_name"] for a in alts),
                "Excluded_Candidates_ClassMismatch": "; ".join(f["resfinder_name"] for f in flagged),
            }
        )

    merged = pd.DataFrame(rows).merge(prev_summary, left_on="ARO_Accession", right_on="ARO Accession", how="left")
    merged = merged.drop(columns=["ARO Accession"])
    merged["max_ncbi_wgs_pct"] = merged["max_ncbi_wgs_pct"].round(3)
    merged = merged.rename(
        columns={
            "n_pathogens_detected_in": "N_Pathogens_With_Hit",
            "max_ncbi_wgs_pct": "Max_NCBI_WGS_Prevalence_Pct",
            "top_pathogen_by_wgs_pct": "Top_Pathogen_By_WGS_Pct",
            "has_perfect_strict_hit": "Has_Perfect_Strict_Hit",
        }
    )
    return merged.sort_values(["Match_Status", "ARO_Accession"]).reset_index(drop=True)


def _json_safe(value):
    """Coerce numpy/pandas scalar types to plain Python types for json.dump,
    and NaN/NaT to None."""
    if pd.isna(value):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def build_json_dictionary(merged: pd.DataFrame) -> dict:
    """Build the JSON-serializable merged dictionary (forward + reverse index).

    Args:
        merged: Output of :func:`build_merged_table`.

    Returns:
        A dict with keys ``metadata``, ``card_aro_to_resfinder`` (CARD ARO
        Accession -> details), and ``resfinder_to_card_aro`` (ResFinder gene
        name -> list of CARD ARO Accessions, matched entries only).
    """
    card_to_resfinder = {}
    for _, r in merged.iterrows():
        card_to_resfinder[r["ARO_Accession"]] = {
            "card_name": _json_safe(r["CARD_Name"]),
            "card_model_type": _json_safe(r["Model_Type"]),
            "card_aro_categories": _json_safe(r["ARO_Categories"]),
            "match_status": _json_safe(r["Match_Status"]),
            "resfinder_gene_name": _json_safe(r["ResFinder_Gene_Name"]) or None,
            "resfinder_class": _json_safe(r["ResFinder_Class"]) or None,
            "match_tier": _json_safe(r["Match_Tier"]) or None,
            "match_confidence": _json_safe(r["Match_Confidence"]) or None,
            "alternate_resfinder_synonyms": [s.strip() for s in r["Alternate_ResFinder_Synonyms"].split(";") if s.strip()],
            "excluded_candidates_class_mismatch": [s.strip() for s in r["Excluded_Candidates_ClassMismatch"].split(";") if s.strip()],
            "n_pathogens_with_hit_ncbi": _json_safe(r["N_Pathogens_With_Hit"]),
            "max_ncbi_wgs_prevalence_pct": _json_safe(r["Max_NCBI_WGS_Prevalence_Pct"]),
            "top_pathogen_by_wgs_pct": _json_safe(r["Top_Pathogen_By_WGS_Pct"]),
        }

    resfinder_to_card: dict[str, list[str]] = {}
    for _, r in merged[merged["Match_Status"] == "matched"].iterrows():
        resfinder_to_card.setdefault(r["ResFinder_Gene_Name"], []).append(r["ARO_Accession"])

    return {
        "metadata": {
            "description": "CARD (RGI wildcard prevalence) ARO-term <-> ResFinder DB gene-name reconciliation",
            "script_version": __version__,
            "n_card_aro_terms_total": len(merged),
            "n_matched": int((merged["Match_Status"] == "matched").sum()),
            "n_no_reliable_match": int((merged["Match_Status"] == "no_reliable_match").sum()),
            "n_out_of_scope_chromosomal_or_regulatory": int(merged["Match_Status"].str.startswith("out_of_scope").sum()),
            "match_tiers_explained": {
                "exact": "identical string",
                "case_insensitive": "identical after case-folding only",
                "loose": "identical after stripping case/spaces/parens/hyphens (apostrophes preserved)",
                "debla": "identical after stripping ResFinder's leading 'bla' beta-lactamase prefix",
                "van_structural": "CARD 'vanX gene in vanY cluster' <-> ResFinder 'VanX-Y' structural rule",
                "alias_bridge": "linked only via an 'Alternative name' cross-reference mined from phenotypes.txt; Medium confidence",
            },
        },
        "card_aro_to_resfinder": card_to_resfinder,
        "resfinder_to_card_aro": resfinder_to_card,
    }


def build_ambiguous_review(match_df: pd.DataFrame) -> pd.DataFrame:
    """Build the transparency table of every CARD entry with >1 raw candidate.

    Args:
        match_df: Output of :func:`apply_manual_overrides`.

    Returns:
        DataFrame for manual audit, one row per multi-candidate CARD entry.
    """
    amb_rows = []
    for _, r in match_df[match_df["n_candidates"] > 1].iterrows():
        prim = r["primary_match"]
        amb_rows.append(
            {
                "ARO_Accession": r["ARO_Accession"],
                "CARD_Name": r["CARD_Name"],
                "ARO_Categories": r["ARO_Categories"],
                "Primary_ResFinder_Match": prim["resfinder_name"] if prim else "(none - all candidates class-inconsistent)",
                "Primary_Tier": prim["tier_tag"] if prim else "",
                "Documented_Synonym_Alternates": "; ".join(
                    f"{m['resfinder_name']} [{m['tier_tag']}]" for m in r["alternate_matches"]
                ),
                "Excluded_ClassMismatch_Candidates": "; ".join(
                    f"{m['resfinder_name']} [{m['tier_tag']}]" for m in r["category_flagged"]
                ),
            }
        )
    return pd.DataFrame(amb_rows).sort_values("ARO_Accession") if amb_rows else pd.DataFrame(
        columns=[
            "ARO_Accession", "CARD_Name", "ARO_Categories", "Primary_ResFinder_Match",
            "Primary_Tier", "Documented_Synonym_Alternates", "Excluded_ClassMismatch_Candidates",
        ]
    )


# ============================================================================
# Command-line interface
# ============================================================================

_EPILOG = """\
usage examples:
  # Defaults (reads from /mnt/user-data/uploads, writes to /mnt/user-data/outputs)
  python3 build_mapping.py

  # Explicit directories
  python3 build_mapping.py --inputdir ./data --outputdir ./results

  # Point at individually renamed/relocated files
  python3 build_mapping.py \\
      --card-prevalence-file /data/card_v3.9/prevalence.tsv \\
      --resfinder-notes-file /data/resfinder_db/notes.txt \\
      --resfinder-phenotypes-file /data/resfinder_db/phenotypes.txt \\
      --outputdir ./results

  # Inspect the effect of the safety-net rules by disabling them
  python3 build_mapping.py --no-category-check --no-manual-overrides -v

Full methodology, validation results, and reference list: see METHODOLOGY.md.
"""


def build_arg_parser() -> argparse.ArgumentParser:
    """Construct the command-line argument parser.

    Returns:
        A fully configured :class:`argparse.ArgumentParser`.
    """
    parser = argparse.ArgumentParser(
        prog="build_mapping.py",
        description=(
            "Build a merged dictionary matching CARD Antibiotic Resistance Ontology "
            "(ARO) terms to ResFinder DB acquired-resistance-gene names, using "
            "literature-grounded deterministic rules (bla-prefix stripping, van-cluster "
            "structural matching, and phenotypes.txt-derived historical-name alias "
            "bridging), plus a data-driven drug-class consistency check that catches "
            "same-spelling/different-gene false positives."
        ),
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )

    io_group = parser.add_argument_group(
        "input / output paths",
        "Any --*-file argument left unset defaults to <inputdir>/<standard filename>.",
    )
    io_group.add_argument(
        "--inputdir",
        default="/mnt/user-data/uploads",
        metavar="DIR",
        help=(
            "Base directory to look for input files in when individual --*-file "
            "arguments are not given. (default: %(default)s)"
        ),
    )
    io_group.add_argument(
        "--outputdir",
        default="/mnt/user-data/outputs",
        metavar="DIR",
        help="Directory to write all output files to; created if it doesn't exist. (default: %(default)s)",
    )
    io_group.add_argument(
        "--card-prevalence-file",
        dest="card_prevalence_file",
        default=None,
        metavar="PATH",
        help=(
            "Path to CARD's RGI 'wildcard' prevalence export (tab-separated; ARO "
            f"Accession, Name, Model Type, Pathogen, NCBI WGS %%, Criteria, ARO "
            f"Categories, ...). (default: <inputdir>/{DEFAULT_FILENAMES['card_prevalence']})"
        ),
    )
    io_group.add_argument(
        "--resfinder-notes-file",
        dest="resfinder_notes_file",
        default=None,
        metavar="PATH",
        help=(
            "Path to ResFinder DB's notes.txt (canonical acquired-gene names grouped "
            f"under '#<Class>:' headers). (default: <inputdir>/{DEFAULT_FILENAMES['resfinder_notes']})"
        ),
    )
    io_group.add_argument(
        "--resfinder-phenotypes-file",
        dest="resfinder_phenotypes_file",
        default=None,
        metavar="PATH",
        help=(
            "Path to ResFinder DB's phenotypes.txt (per-reference-sequence phenotype "
            "table; its free-text Notes column is mined for 'Alternative name' "
            f"cross-references). (default: <inputdir>/{DEFAULT_FILENAMES['resfinder_phenotypes']})"
        ),
    )
    io_group.add_argument(
        "--output-prefix",
        default="card_resfinder_merged_dictionary",
        metavar="NAME",
        help="Base filename (without extension) for the primary CSV/JSON outputs. (default: %(default)s)",
    )

    behavior_group = parser.add_argument_group("pipeline behavior")
    behavior_group.add_argument(
        "--no-category-check",
        action="store_true",
        help=(
            "Disable the data-driven drug-class consistency re-ranking step. Without "
            "it, the highest raw tier always wins, even when it's a known false "
            "positive like CARD 'FAR-1' (a beta-lactamase) string-matching ResFinder's "
            "unrelated 'far1' (fusidic-acid resistance). Useful for auditing what the "
            "safety net actually changes."
        ),
    )
    behavior_group.add_argument(
        "--no-manual-overrides",
        action="store_true",
        help=(
            "Disable the small table of literature-grounded manual curation overrides "
            "(currently one entry: APH(3'')-Ib -> strA). See MANUAL_OVERRIDES in this "
            "file for the documented rationale."
        ),
    )
    behavior_group.add_argument(
        "--save-intermediate",
        action="store_true",
        help=(
            "Also save intermediate pandas objects (as pickles, under "
            "<outputdir>/intermediate/) for reuse in a follow-up Python session "
            "without re-running the full pipeline."
        ),
    )

    verbosity_group = parser.add_argument_group("console output").add_mutually_exclusive_group()
    verbosity_group.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print per-entry diagnostic detail (every manual override and every drug-class-mismatch flag), in addition to the summary statistics.",
    )
    verbosity_group.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress all non-error console output.",
    )

    return parser


def resolve_input_paths(args: argparse.Namespace) -> argparse.Namespace:
    """Fill in default file paths (relative to --inputdir) for any --*-file
    argument the user didn't explicitly set.

    Args:
        args: Parsed arguments from :func:`build_arg_parser`.

    Returns:
        ``args``, mutated in place, with all three ``*_file`` attributes
        guaranteed to be non-None.
    """
    if args.card_prevalence_file is None:
        args.card_prevalence_file = os.path.join(args.inputdir, DEFAULT_FILENAMES["card_prevalence"])
    if args.resfinder_notes_file is None:
        args.resfinder_notes_file = os.path.join(args.inputdir, DEFAULT_FILENAMES["resfinder_notes"])
    if args.resfinder_phenotypes_file is None:
        args.resfinder_phenotypes_file = os.path.join(args.inputdir, DEFAULT_FILENAMES["resfinder_phenotypes"])
    return args


def _validate_input_files(args: argparse.Namespace) -> None:
    """Exit with a clear message if any required input file is missing."""
    for label, path in [
        ("--card-prevalence-file", args.card_prevalence_file),
        ("--resfinder-notes-file", args.resfinder_notes_file),
        ("--resfinder-phenotypes-file", args.resfinder_phenotypes_file),
    ]:
        if not os.path.isfile(path):
            sys.exit(
                f"ERROR: {label} not found: '{path}'\n"
                f"       Set --inputdir to the directory containing your CARD/ResFinder "
                f"files, or pass {label} explicitly."
            )


def main(argv: list[str] | None = None) -> None:
    """Entry point: parse arguments, run the pipeline, write outputs.

    Args:
        argv: Argument list to parse instead of ``sys.argv[1:]`` (mainly for
            testing); ``None`` uses the real command-line arguments.
    """
    args = build_arg_parser().parse_args(argv)
    args = resolve_input_paths(args)
    _validate_input_files(args)
    os.makedirs(args.outputdir, exist_ok=True)

    def log(msg: str) -> None:
        if not args.quiet:
            print(msg)

    def log_verbose(msg: str) -> None:
        if args.verbose and not args.quiet:
            print(msg)

    # ---- 1. Load ----------------------------------------------------------
    log(f"Loading CARD prevalence data from {args.card_prevalence_file} ...")
    card, card_unique = load_card_prevalence(args.card_prevalence_file)
    log(f"  {len(card):,} rows -> {len(card_unique):,} unique ARO terms")

    log(f"Loading ResFinder notes from {args.resfinder_notes_file} ...")
    notes_df = parse_resfinder_notes(args.resfinder_notes_file)
    log(f"  {len(notes_df):,} unique ResFinder gene names across {notes_df['notes_class'].nunique()} classes")

    log(f"Loading ResFinder phenotypes from {args.resfinder_phenotypes_file} ...")
    pheno = load_resfinder_phenotypes(args.resfinder_phenotypes_file)
    n_alt = sum(len(a) for a in pheno["alt_names"])
    log(f"  {len(pheno):,} rows; {n_alt} 'Alternative name' cross-references mined")

    # ---- 2. Alias graph + normalization index ------------------------------
    alias_uf = build_alias_graph(pheno)
    resfinder_df = build_resfinder_alias_table(notes_df, alias_uf)
    rf_index = build_resfinder_lookup_index(resfinder_df)

    # ---- 3. Match -----------------------------------------------------------
    log("Matching CARD ARO terms to ResFinder gene names ...")
    match_df = match_card_to_resfinder(card_unique, resfinder_df, rf_index)

    match_df = apply_category_consistency_check(match_df, enabled=not args.no_category_check)
    if not args.no_category_check:
        for _, r in match_df.iterrows():
            if r["category_flagged"]:
                flagged_names = [m["resfinder_name"] for m in r["category_flagged"]]
                prim = r["primary_match"]["resfinder_name"] if r["primary_match"] else None
                log_verbose(
                    f"  class-mismatch flagged: {r['ARO_Accession']} '{r['CARD_Name']}' "
                    f"[{r['ARO_Categories']}] primary={prim} excluded={flagged_names}"
                )

    match_df = apply_manual_overrides(
        match_df, MANUAL_OVERRIDES, enabled=not args.no_manual_overrides, log=log_verbose
    )

    # ---- 4. Summary stats -----------------------------------------------------
    has_primary = match_df["primary_match"].apply(lambda m: m is not None)
    log(f"\nCARD unique ARO terms total: {len(card_unique)}")
    log(f"  -> matched to a confident ResFinder gene: {int(has_primary.sum())}")
    log(f"  -> no reliable ResFinder match / out of scope: {int((~has_primary).sum())}")
    log(f"  -> had >1 raw candidate before class-consistency rerank: {int((match_df['n_candidates'] > 1).sum())}")

    tier_counts: Counter = Counter()
    for _, r in match_df.iterrows():
        for m in r["matches"]:
            tier_counts[m["tier_tag"]] += 1
    log("\nMatch tier distribution (counts across all match links):")
    for tier in TIER_PRIORITY:
        log(f"  {tier}: {tier_counts.get(tier, 0)}")

    # ---- 5. Prevalence context + final tables ----------------------------------
    prev_summary = summarize_prevalence(card)
    merged = build_merged_table(match_df, prev_summary)
    output_json = build_json_dictionary(merged)
    ambiguous_df = build_ambiguous_review(match_df)

    # ---- 6. Write outputs -----------------------------------------------------
    csv_path = os.path.join(args.outputdir, f"{args.output_prefix}.csv")
    json_path = os.path.join(args.outputdir, f"{args.output_prefix}.json")
    ambiguous_path = os.path.join(args.outputdir, "ambiguous_cases_review.csv")

    merged.to_csv(csv_path, index=False)
    with open(json_path, "w") as fh:
        json.dump(output_json, fh, indent=2)
    ambiguous_df.to_csv(ambiguous_path, index=False)

    if args.save_intermediate:
        inter_dir = os.path.join(args.outputdir, "intermediate")
        os.makedirs(inter_dir, exist_ok=True)
        match_df.to_pickle(os.path.join(inter_dir, "match_df.pkl"))
        resfinder_df.to_pickle(os.path.join(inter_dir, "resfinder_df.pkl"))
        card_unique.to_pickle(os.path.join(inter_dir, "card_unique.pkl"))
        log(f"\nIntermediate pickles saved to {inter_dir}/")

    log(f"\nOutputs written to {args.outputdir}/:")
    log(f"  {os.path.basename(csv_path)}")
    log(f"  {os.path.basename(json_path)}")
    log(f"  {os.path.basename(ambiguous_path)}")


if __name__ == "__main__":
    main()
