"""
Data loading and preparation for the Crystal Explorer viewer.

Reads Grouped_conditions.csv -- the file plot.py (the pipeline's stage-3
script) produces per protein, alongside the static PNG plots and PDF
table. That file is *already* grouped (one row per unique condition, PDB
IDs merged) and already carries a `row_id` column, so this module's job
is just:

  1. Build an interactive Plotly scatter (pH vs PEG concentration), each
     point tagged with its `row_id` in `customdata`.
  2. Build the matching HTML table rows, each carrying the same `row_id`
     as its DOM id `row-<row_id>`.

Using the same Grouped_conditions.csv that produced the static PDF/PNG
keeps the web app and the static exports showing identical groupings and
identical row numbering.

Column expectations (from plot.py's run_plot -> Grouped_conditions.csv):
    row_id, PDB_ID, Score, Seq_id, Pubmed_id, Polymer, Assembly, Method,
    Ligands, pH, Temp, plot_pH_numeric, PEG_con_plot, compound

Missing columns are tolerated -- anything not present is filled with NaN
so the page still renders instead of erroring out.
"""

import ast
import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go

REQUIRED_COLUMNS = [
    "row_id", "PDB_ID", "UniProt", "UniProt_query", "Score", "Seq_id", "Pubmed_id", "Polymer", "Assembly",
    "Method", "Non_polymers", "plot_pH_numeric", "Temp", "PEG_con_plot", "compound",
]

# Method -> Plotly marker symbol (mirrors plot.py's matplotlib marker mapping)
METHOD_MARKER_MAP = {
    "vapor diffusion hanging drop": "square",
    "hanging drop": "square",
    "hanging drop vapor diffusion": "square",
    "vapor diffusion sitting drop": "triangle-up",
    "sitting drop": "triangle-up",
    "sitting drop vapor diffusion": "triangle-up",
    "vapor diffusion": "circle",
    "unspecified": "x",
    "": "x",
}
FALLBACK_MARKERS = ["diamond", "triangle-down", "pentagon", "star", "triangle-left", "triangle-right"]

NO_PH_Y = 4.0         # "No pH" row position
NO_PEG_X = -5.0       # "No PEG" column position
NO_TEMP_X = 270.0     # "No Temp" column position (below the ~275K typical minimum)

# Score color scale shared by both plots' markers and the HTML legend below
# them -- kept as plain hex stops (rather than relying on Plotly's internal
# "Viridis" table) so the CSS gradient in the legend matches the marker
# colors exactly.
SCORE_MIN = 0.5
SCORE_MAX = 1.0
SCORE_TICKS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
VIRIDIS_STOPS = [
    "#440154", "#482878", "#3e4989", "#31688e", "#26828e",
    "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725",
]


def score_gradient_stops() -> str:
    """Color stops (no direction) matching the marker colorscale (cmin=0.5,
    cmax=1.0) -- fed into a CSS custom property so the template can point
    the gradient horizontally or vertically as the layout needs."""
    n = len(VIRIDIS_STOPS)
    return ", ".join(f"{color} {i / (n - 1) * 100:.1f}%" for i, color in enumerate(VIRIDIS_STOPS))


def _ensure_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan
    return df


def _assign_marker(method: str, used_fallback: dict) -> str:
    method = str(method).lower().strip()
    if method in METHOD_MARKER_MAP:
        return METHOD_MARKER_MAP[method]
    if method not in used_fallback:
        used_fallback[method] = len(used_fallback) % len(FALLBACK_MARKERS)
    return FALLBACK_MARKERS[used_fallback[method]]


def _method_symbol_map(df: pd.DataFrame) -> dict:
    """method -> marker symbol, assigned once (in sorted method order) so
    the PEG plot, Temp plot, and the shared HTML legend all agree -- this
    is what lets both plots share one "Method" legend."""
    used_fallback: dict = {}
    return {method: _assign_marker(method, used_fallback) for method in sorted(df["Method"].unique(), key=str)}


def build_method_legend(df: pd.DataFrame):
    """List of {"method": <title-cased label>, "symbol": <plotly marker
    symbol>} covering every method in df -- shared by both plots since
    they're built from the same df."""
    entries = []
    for method, symbol in _method_symbol_map(df).items():
        label = str(method).title() if str(method).strip() else "Unspecified"
        entries.append({"method": label, "symbol": symbol})
    return entries


def load_conditions(grouped_csv_path: str) -> pd.DataFrame:
    """Load Grouped_conditions.csv (already grouped + row_id'd by plot.py)
    and coerce numeric columns -- no re-grouping needed here."""
    df = pd.read_csv(grouped_csv_path)
    df = _ensure_columns(df)

    df["Method"] = df["Method"].fillna("unspecified").astype(str)
    df.loc[df["Method"] == "", "Method"] = "unspecified"

    for col in ["Score", "Seq_id", "Temp", "plot_pH_numeric", "PEG_con_plot"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["Score"] = df["Score"].clip(SCORE_MIN, SCORE_MAX)

    df = df.sort_values("Score", ascending=False).reset_index(drop=True)
    return df


def _build_scatter_plot(df: pd.DataFrame, x_col: str, x_fallback: float,
                         x_label: str, x_hover_label: str, x_hover_fmt: str) -> go.Figure:
    """Generic pH-vs-<x_col> scatter plot builder, shared by the PEG and
    Temperature plots. Each point's customdata is its row_id, used by the
    front-end click handler to scroll to / highlight the matching table row."""
    fig = go.Figure()

    if df.empty:
        fig.update_layout(title="No data to display", xaxis_title=x_label, yaxis_title="pH")
        return fig

    symbol_map = _method_symbol_map(df)
    for method, group in df.groupby("Method"):
        marker_symbol = symbol_map[method]

        x = group[x_col].fillna(x_fallback)
        y = group["plot_pH_numeric"].fillna(NO_PH_Y)
        score = group["Score"].fillna(0.5)

        hover_text = []
        for rec in group.to_dict("records"):
            score_str = f"{rec['Score']:.3f}" if pd.notna(rec["Score"]) else "n/a"
            ph_str = f"{rec['plot_pH_numeric']:.2f}" if pd.notna(rec["plot_pH_numeric"]) else "n/a"
            x_str = x_hover_fmt.format(rec[x_col]) if pd.notna(rec[x_col]) else "n/a"
            hover_text.append(
                f"<b>{rec['PDB_ID']}</b><br>"
                f"Method: {rec['Method']}<br>"
                f"pH: {ph_str}<br>"
                f"{x_hover_label}: {x_str}<br>"
                f"Score: {score_str}"
            )

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="markers",
            name=method.title() if method else "Unspecified",
            showlegend=False,
            marker=dict(
                symbol=marker_symbol,
                size=14,
                color=score,
                colorscale="Viridis",
                cmin=SCORE_MIN,
                cmax=SCORE_MAX,
                line=dict(width=1, color="#14181B"),
                showscale=False,
            ),
            customdata=group["row_id"],
            text=group["PDB_ID"],
            hovertext=hover_text,
            hoverinfo="text",
        ))

    # No per-plot legend or colorbar -- both plots share one "Method"
    # legend and one Score gradient, rendered as HTML below the plots
    # (see build_method_legend / score_gradient_stops), so each plot's own
    # space stays a plain square instead of losing width to a side legend.
    fig.update_layout(
        xaxis_title=x_label,
        yaxis_title="pH",
        showlegend=False,
        template="plotly_white",
        font=dict(family="IBM Plex Sans, sans-serif", size=13, color="#14181B"),
        margin=dict(l=60, r=20, t=20, b=50),
        autosize=True,
        clickmode="event+select",
    )
    fig.update_xaxes(zeroline=False, gridcolor="#E1E5E8")
    fig.update_yaxes(zeroline=False, gridcolor="#E1E5E8")

    return fig


def build_peg_plot(df: pd.DataFrame) -> go.Figure:
    """pH vs PEG concentration (%)."""
    return _build_scatter_plot(df, "PEG_con_plot", NO_PEG_X, "PEG concentration (%)",
                                "PEG", "{:.1f}%")


def build_temp_plot(df: pd.DataFrame) -> go.Figure:
    """pH vs Temperature (K)."""
    return _build_scatter_plot(df, "Temp", NO_TEMP_X, "Temperature (K)",
                                "Temp", "{:.0f}K")


# Buffer / Salt / Precipitant / Additive -- kept in this fixed order
# everywhere (legend, subplot grid) so the four categories always land in
# the same place regardless of which one happens to be biggest that run.
CATEGORY_ORDER = ["Buffer", "Salt", "Precipitant", "Additive"]
CATEGORY_COLORS = {
    "Buffer": "#0B6E7D",       # matches the app's --accent teal
    "Salt": "#E69F00",
    "Precipitant": "#0072B2",
    "Additive": "#8A94A6",
}


def build_reagent_frequency(df: pd.DataFrame) -> list:
    """Aggregate reagent frequency and typical concentration across every
    row's parsed 'compound' cell.

    Returns a list of dicts, one per distinct reagent -- {"name",
    "category", "count", "total", "top_concentration",
    "top_concentration_count"} -- sorted by count descending. `count` is
    the number of *conditions* (rows) the reagent appears in, not total
    mentions, so a reagent occurring twice within one row's clause list
    (rare, but possible) still counts once. `top_concentration` is the
    single most common exact concentration string seen for that reagent
    (the mode) rather than a numeric average, since different rows can
    report the same reagent in different units (%, M, mM) that can't be
    meaningfully averaged together -- the mode stays an exact, directly
    reproducible value instead."""
    from collections import Counter
    from pipeline.compound_extraction import get_reagent_category

    total = len(df)
    hit_counts = Counter()
    concentration_counts = {}

    for cell in df.get("compound", pd.Series(dtype=str)):
        seen_in_row = set()
        for entry in _parse_compound_cell(cell):
            name = entry["name"]
            if name not in seen_in_row:
                hit_counts[name] += 1
                seen_in_row.add(name)
            concentration_counts.setdefault(name, Counter())[entry["concentration"]] += 1

    results = []
    for name, count in hit_counts.items():
        top_conc, top_conc_count = concentration_counts[name].most_common(1)[0]
        results.append({
            "name": name,
            "category": get_reagent_category(name),
            "count": count,
            "total": total,
            "top_concentration": top_conc,
            "top_concentration_count": top_conc_count,
        })
    results.sort(key=lambda r: r["count"], reverse=True)
    return results


def _reagent_hover_text(entries: list) -> list:
    return [
        f"<b>{e['name']}</b><br>"
        f"{e['count']}/{e['total']} conditions ({100 * e['count'] / e['total']:.0f}%)<br>"
        f"Most common: {e['top_concentration']} "
        f"({e['top_concentration_count']} of {e['count']})"
        for e in entries
    ]


def build_reagent_frequency_plot(freq: list, top_n: int = 15) -> go.Figure:
    """Flat, ranked-by-frequency horizontal bar chart across all reagents
    (top `top_n`), colored by category so the buffer/salt/precipitant/
    additive split is still visible even in the combined view."""
    fig = go.Figure()

    if not freq:
        fig.update_layout(title="No compounds extracted", template="plotly_white")
        return fig

    top = freq[:top_n]
    # Plotly lists horizontal bars bottom-to-top, so reverse to put the
    # highest count at the top of the chart, matching reading order.
    top = list(reversed(top))

    for category in CATEGORY_ORDER:
        entries = [e for e in top if e["category"] == category]
        if not entries:
            continue
        fig.add_trace(go.Bar(
            x=[e["count"] for e in entries],
            y=[e["name"] for e in entries],
            orientation="h",
            name=category,
            marker=dict(color=CATEGORY_COLORS[category]),
            hovertext=_reagent_hover_text(entries),
            hoverinfo="text",
        ))

    fig.update_layout(
        barmode="overlay",
        xaxis_title="Conditions",
        template="plotly_white",
        font=dict(family="IBM Plex Sans, sans-serif", size=13, color="#14181B"),
        margin=dict(l=180, r=20, t=20, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        autosize=True,
    )
    fig.update_xaxes(zeroline=False, gridcolor="#E1E5E8")
    fig.update_yaxes(zeroline=False)
    return fig


def build_reagent_frequency_by_category_plot(freq: list, top_n: int = 8) -> go.Figure:
    """2x2 grid -- one horizontal bar chart per Buffer/Salt/Precipitant/
    Additive category, each showing its own top `top_n` reagents. Mirrors
    how a crystallization cocktail actually gets assembled (one buffer +
    one salt + one precipitant, plus optional additives), so this is the
    more directly "what should I pick" view versus the flat overall one."""
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=2, cols=2, subplot_titles=CATEGORY_ORDER,
        horizontal_spacing=0.15, vertical_spacing=0.18,
    )
    positions = {"Buffer": (1, 1), "Salt": (1, 2), "Precipitant": (2, 1), "Additive": (2, 2)}

    has_any = False
    for category in CATEGORY_ORDER:
        entries = [e for e in freq if e["category"] == category][:top_n]
        if not entries:
            continue
        has_any = True
        entries = list(reversed(entries))
        row, col = positions[category]
        fig.add_trace(go.Bar(
            x=[e["count"] for e in entries],
            y=[e["name"] for e in entries],
            orientation="h",
            marker=dict(color=CATEGORY_COLORS[category]),
            hovertext=_reagent_hover_text(entries),
            hoverinfo="text",
            showlegend=False,
        ), row=row, col=col)

    if not has_any:
        fig.update_layout(title="No compounds extracted", template="plotly_white")
        return fig

    fig.update_layout(
        template="plotly_white",
        font=dict(family="IBM Plex Sans, sans-serif", size=13, color="#14181B"),
        margin=dict(l=20, r=20, t=40, b=20),
        autosize=True,
        height=560,
    )
    fig.update_xaxes(zeroline=False, gridcolor="#E1E5E8", title_text="Conditions")
    fig.update_yaxes(zeroline=False)
    return fig


# Kept for backward compatibility with any existing callers.
build_plot = build_peg_plot


TABLE_COLUMNS = [
    ("PDB_ID", "PDB ID"),
    ("UniProt", "UniProt"),
    ("Score", "Score"),
    ("Seq_id", "Seq. ID (%)"),
    ("Pubmed_id", "PubMed"),
    ("Polymer", "Polymer"),
    ("Assembly", "Assembly"),
    ("Method", "Method"),
    ("View3D", "3D View"),
    ("pH", "pH"),
    ("Temp", "Temp (K)"),
    ("compound", "Compounds"),
]


def _format_pubmed(value) -> str:
    # RCSB's own API returns pdbx_database_id_pub_med as -1 (not null/"NA")
    # for entries that simply have no PubMed ID -- treat that the same as
    # missing rather than showing/linking a bogus "-1".
    if pd.isna(value) or str(value) in ("NA", "", "-1", "-1.0"):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _parse_non_polymer_codes(value) -> str:
    """Non_polymers cells hold a dict -- {"non_polymers": "LIG1, LIG2",
    "view_3d_url": "https://www.rcsb.org/3d-view/1XOM"} -- that comes back
    as a stringified dict after the CSV round-trip. Pull out just the
    ligand codes, if any; "" if there's nothing usable (blank cell,
    unparsable)."""
    if pd.isna(value) or not str(value).strip():
        return ""
    parsed = value if isinstance(value, dict) else None
    if parsed is None:
        try:
            parsed = ast.literal_eval(str(value))
        except (ValueError, SyntaxError):
            return ""
    if not isinstance(parsed, dict):
        return ""
    return parsed.get("non_polymers") or ""


_COMPOUND_CELL_RE = re.compile(r"'([^']+)'\s*\(([^)]*)\)")


def _parse_compound_cell(value) -> list:
    """The 'compound' cell holds compound_extraction.py's
    format_compounds() output -- "'PEG 4000' (25.0 % w/v), 'NaCl' (0.1 M)"
    -- one or more single-quoted compound names each followed by a
    parenthesized concentration. Split back into (name, concentration)
    pairs so each name can link to its own PubChem page instead of the
    whole cell being unclickable free text."""
    if pd.isna(value) or not str(value).strip():
        return []
    return [
        {"name": name.strip(), "concentration": concentration.strip()}
        for name, concentration in _COMPOUND_CELL_RE.findall(str(value))
    ]


def _split_merged_ids(value) -> list:
    """PDB_ID and UniProt cells hold plot.py's merged ", "-joined string
    when several structures share one crystallization condition (e.g.
    "1ABC, 2XYZ"). Split back into individual IDs so each can get its own
    link -- a single link wrapping the whole merged string is not a valid
    PDB ID/accession and can't be clicked "separately" for each entry."""
    if pd.isna(value) or not str(value).strip():
        return []
    return [pid.strip() for pid in str(value).split(",") if pid.strip()]


# Small colorblind-safe categorical palette (Okabe-Ito), cycled by index.
# Only ever used to link a specific merged PDB_ID to its own UniProt
# accession(s) within one row -- meaningless across rows, so 6 colors is
# plenty even though a row could in principle merge more PDB_IDs than that.
_UNIPROT_COLORS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2"]


def _split_per_pdb(value) -> list:
    """UniProt/UniProt_query cells hold plot.py's " | "-joined string --
    one segment per merged PDB_ID, in the same order as the PDB_ID cell,
    each segment itself a ", "-joined list of accessions for that one
    PDB_ID (e.g. "P0DTD1, P05161 | P0DTD1, J3QS39" for two merged
    complexes). Returns a list of per-PDB_ID accession lists. A row from
    before this format existed has no "|" at all -- comes back as a
    single segment, which build_table_rows() treats as unaligned/legacy."""
    if pd.isna(value) or not str(value).strip():
        return []
    return [_split_merged_ids(segment) for segment in str(value).split("|")]


def _build_3d_view(pdb_ids, non_polymers_value):
    """Every row has at least one PDB ID, so every row gets 3D-view link(s)
    -- unlike the old ligand-only link, this doesn't depend on the entry
    having non-polymers. Ligand codes (if any) are still surfaced as the
    link label so that information isn't lost.

    A row can carry several merged PDB IDs (see `_split_merged_ids`), and each
    structure needs its own link/box -- one link for the whole group would
    point at a single bogus molstar URL and couldn't open a specific
    structure. So this returns one dict per PDB ID.

    Each link points at molstar.org's bare viewer (?pdb=<id>) rather than
    RCSB's own https://www.rcsb.org/3d-view/<id> page: the RCSB page is a
    full site (header, tabs, summary text) meant to be browsed on its own,
    not embedded, whereas molstar.org/viewer is built to be embedded --
    just the 3D canvas and its control panel, which is what the in-page
    viewer panel here wants to show."""
    codes = _parse_non_polymer_codes(non_polymers_value)
    return [
        {
            "url": f"https://molstar.org/viewer/?pdb={pdb_id}",
            "pdb_id": pdb_id,
            "codes": codes,
        }
        for pdb_id in pdb_ids
    ]


def build_table_rows(df: pd.DataFrame):
    """Return a list of dicts: {"row_id": ..., "cells": [(key, value), ...]}
    with cells in the same order as TABLE_COLUMNS, ready for simple
    template iteration. Cells are (key, value) pairs -- rather than bare
    values -- so the template can special-case columns like "View3D"
    (rendered as a clickable link) without extra template logic elsewhere."""
    rows = []
    for rec in df.to_dict("records"):
        pdb_ids = _split_merged_ids(rec["PDB_ID"])
        uniprot_per_pdb = _split_per_pdb(rec["UniProt"])
        # UniProt_query is also " | "-joined per PDB_ID now (see plot.py) --
        # flattened into one set since a match against *any* merged PDB
        # entry's query-matched accession is what "is_query_match" means.
        query_uniprot_ids = {acc for segment in _split_per_pdb(rec["UniProt_query"]) for acc in segment}

        # Only color-code PDB_ID <-> UniProt when there's real ambiguity to
        # resolve: more than one merged PDB_ID, the per-PDB_ID accession
        # data actually lines up with them (older rows predate the " | "
        # format and come back as one unaligned segment), and the
        # accessions genuinely differ between at least two of them -- a
        # uniform row (every PDB_ID -> the same single UniProt) has
        # nothing to disambiguate, so it renders exactly as before.
        aligned = len(uniprot_per_pdb) == len(pdb_ids)
        distinct_sets = {tuple(seg) for seg in uniprot_per_pdb} if aligned else set()
        use_colors = len(pdb_ids) > 1 and aligned and len(distinct_sets) > 1

        if use_colors:
            pdb_id_cells = [
                {"id": pid, "color": _UNIPROT_COLORS[i % len(_UNIPROT_COLORS)]}
                for i, pid in enumerate(pdb_ids)
            ]
            uniprot_groups = [
                {
                    "color": _UNIPROT_COLORS[i % len(_UNIPROT_COLORS)],
                    "accessions": [
                        {"accession": acc, "is_query_match": acc in query_uniprot_ids}
                        for acc in uniprot_per_pdb[i]
                    ],
                }
                for i in range(len(pdb_ids))
                if uniprot_per_pdb[i]
            ]
        else:
            pdb_id_cells = [{"id": pid, "color": None} for pid in pdb_ids]
            seen = []
            for segment in uniprot_per_pdb:
                for acc in segment:
                    if acc not in seen:
                        seen.append(acc)
            uniprot_groups = [{
                "color": None,
                "accessions": [
                    {"accession": acc, "is_query_match": acc in query_uniprot_ids}
                    for acc in seen
                ],
            }] if seen else []

        values = {
            "PDB_ID": pdb_id_cells,
            "UniProt": uniprot_groups,
            "Score": f"{rec['Score']:.3f}" if pd.notna(rec["Score"]) else "",
            "Seq_id": f"{rec['Seq_id']:.1f}" if pd.notna(rec["Seq_id"]) else "",
            "Pubmed_id": _format_pubmed(rec["Pubmed_id"]),
            "Polymer": rec["Polymer"] if pd.notna(rec["Polymer"]) else "",
            "Assembly": rec["Assembly"] if pd.notna(rec["Assembly"]) else "",
            "Method": rec["Method"].title() if isinstance(rec["Method"], str) else "",
            "View3D": _build_3d_view(pdb_ids, rec["Non_polymers"]),
            "pH": f"{rec['plot_pH_numeric']:.2f}" if pd.notna(rec["plot_pH_numeric"]) else "",
            "Temp": f"{rec['Temp']:.1f}" if pd.notna(rec["Temp"]) else "",
            "compound": _parse_compound_cell(rec["compound"]),
        }
        rows.append({
            "row_id": rec["row_id"],
            "cells": [(key, values[key]) for key, _ in TABLE_COLUMNS],
        })
    return rows