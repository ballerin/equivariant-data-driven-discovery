#!/usr/bin/env python3
"""
Generate LaTeX tables of the pointwise -> weak collinearity penalty.

For every (equation, library, noise level) cell the percent change

    pct = 100 * (weak - pointwise) / pointwise

is computed from a collinearity diagnostics CSV, then averaged two ways:

    Table 1  per equation, across libraries.  The direct and IBP weak
             forms are pooled into a single mean/median, since IBP is a
             library-level choice rather than an equation-level one.
    Table 2  per library, across equations.  The direct and IBP weak
             forms get one row each, so the effect of integration by
             parts is visible library by library.

Both the arithmetic mean and the median are reported.  The distribution
of ``pct`` spans several orders of magnitude, so the median is the more
representative of the two; the mean is retained for completeness.

Usage
-----
    python scripts/build_collinearity_latex.py \\
        results/collinearity_diagnostics_spectral.csv

    # kappa_norm instead of max VIF, written to a file
    python scripts/build_collinearity_latex.py results/coll.csv \\
        --metric cond_normalised -o tables.tex

    # both metrics, dropping rows where max_vif hit the solver cap
    python scripts/build_collinearity_latex.py results/coll.csv \\
        --metric both --drop-censored
"""

import argparse
import csv
import os
import statistics as st
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
sys.path.insert(0, SRC_DIR)

from equivariant_pde_find.summary_metadata import (  # noqa: E402
    EQ_NAMES, normalize_sweep_row,
)

# ── constants ──────────────────────────────────────────────────────────

#: ``max_vif`` is censored at this value by the diagnostics script; rows
#: sitting on it are lower bounds, not measurements.
VIF_CAP = 1e10

WEAK_FORMS = ["space_weak_direct", "space_weak_ibp"]

EQ_ORDER = ["ac", "ch", "ks_cell", "ks_wchaos", "pme", "sh"]
LIB_ORDER = ["multigraph", "mod_multigraph", "pdefind", "wsindy", "polyjet"]
NOISE_ORDER = ["clean", "noise2pct"]

LIB_NAMES = {
    "multigraph": "Multigraphs",
    "mod_multigraph": "Modified multigraphs",
    "pdefind": "PDE-FIND",
    "wsindy": "WSINDy",
    "polyjet": "PolyJet",
}

NOISE_NAMES = {"clean": "Clean", "noise2pct": r"2\% noise"}

METRIC_NAMES = {
    "max_vif": r"$\mathrm{VIF}_{\max}$",
    "mean_vif": r"$\overline{\mathrm{VIF}}$",
    "cond_normalised": r"$\kappa_{\mathrm{norm}}$",
}

NUMERIC_COLS = (
    "max_abs_corr", "max_vif", "mean_vif", "cond_raw", "cond_normalised",
)

# ── data loading ───────────────────────────────────────────────────────


def load_rows(csv_path):
    """Return diagnostics rows keyed by (eq, noise, library, formulation)."""
    with open(csv_path, newline="") as fh:
        rows = [normalize_sweep_row(r) for r in csv.DictReader(fh)]
    for r in rows:
        for k in NUMERIC_COLS:
            if r.get(k) not in (None, ""):
                r[k] = float(r[k])
    return {
        (r["eq_name"], r["noise_level"], r["library"], r["formulation"]): r
        for r in rows
    }


def pct_changes(data, metric, weak_forms, axis, drop_censored=False):
    """Percent changes pointwise -> weak, grouped along ``axis``.

    ``weak_forms`` is a list; passing more than one pools them into the
    same group.  ``axis`` is ``"eq_name"`` or ``"library"``.  Returns a
    dict mapping ``(group, noise_level)`` to a list of percent changes.
    """
    out = {}
    for (eq, noise, lib, form), pw in data.items():
        if form != "pointwise":
            continue
        group = eq if axis == "eq_name" else lib
        for weak_form in weak_forms:
            weak = data.get((eq, noise, lib, weak_form))
            if weak is None:
                continue
            if (drop_censored and metric == "max_vif"
                    and weak[metric] >= VIF_CAP):
                continue
            pct = 100.0 * (weak[metric] - pw[metric]) / pw[metric]
            out.setdefault((group, noise), []).append(pct)
    return out

# ── LaTeX emission ─────────────────────────────────────────────────────


def fmt_pct(value):
    """Format a percent change for LaTeX, with thin-space thousands."""
    if abs(value) < 10:
        body = f"{abs(value):.1f}"
    else:
        body = f"{abs(value):,.0f}".replace(",", r"\,")
    sign = "+" if value >= 0 else "$-$"
    return rf"{sign}{body}\%"


def _row(stats, group, display):
    """Render one table row: mean and median per noise level."""
    cells = [display]
    for noise in NOISE_ORDER:
        pcts = stats.get((group, noise))
        if not pcts:
            cells += ["---"] * 2
            continue
        cells += [fmt_pct(st.mean(pcts)), fmt_pct(st.median(pcts))]
    return " & ".join(cells) + r" \\"


def build_table(data, metric, axis, drop_censored):
    """Return a LaTeX ``table`` float summarising one grouping axis."""
    metric_name = METRIC_NAMES.get(metric, metric)

    if axis == "eq_name":
        # IBP is a library-level choice, so pool both weak forms here.
        stats = pct_changes(
            data, metric, WEAK_FORMS, axis, drop_censored,
        )
        groups = [g for g in EQ_ORDER if any(k[0] == g for k in stats)]
        head, avg_over = "Equation", "libraries"
        caption_extra = (
            r"  Both the direct and integration-by-parts weak forms are "
            r"pooled into each average."
        )
        body = [_row(stats, g, EQ_NAMES.get(g, g)) for g in groups]
    else:
        per_form = {
            wf: pct_changes(data, metric, [wf], axis, drop_censored)
            for wf in WEAK_FORMS
        }
        stats = per_form["space_weak_direct"]
        groups = [g for g in LIB_ORDER if any(k[0] == g for k in stats)]
        head, avg_over = "Library", "equations"
        caption_extra = (
            r"  Each library is reported for the direct weak form and, "
            r"below it, the integration-by-parts variant."
        )
        body = []
        for g in groups:
            name = LIB_NAMES.get(g, g.replace("_", r"\_"))
            body.append(_row(per_form["space_weak_direct"], g, name))
            if any(k[0] == g for k in per_form["space_weak_ibp"]):
                # Indented under its parent, so the name need not repeat.
                body.append(_row(per_form["space_weak_ibp"], g, r"\quad (IBP)"))

    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        rf"\caption{{Percent change in {metric_name} when moving from the "
        rf"pointwise to the weak formulation, averaged over {avg_over}.  "
        rf"Positive values mean the weak formulation is more collinear."
        rf"{caption_extra}}}",
        rf"\label{{tab:coll-pct-"
        rf"{'eq' if axis == 'eq_name' else 'lib'}-"
        rf"{metric.replace('_', '-')}}}",
        r"\begin{tabular}{l" + "rr" * len(NOISE_ORDER) + "}",
        r"\toprule",
    ]

    # grouped header: one block of two columns per noise level
    top, rule, col = [head], [], 2
    for noise in NOISE_ORDER:
        top.append(rf"\multicolumn{{2}}{{c}}{{{NOISE_NAMES[noise]}}}")
        rule.append(rf"\cmidrule(lr){{{col}-{col + 1}}}")
        col += 2
    lines.append(" & ".join(top) + r" \\")
    lines.append("".join(rule))
    lines.append(
        " & ".join([""] + [r"Mean & Median"] * len(NOISE_ORDER)) + r" \\"
    )
    lines.append(r"\midrule")
    lines += body
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


def build_document(csv_path, metrics, drop_censored):
    data = load_rows(csv_path)
    out = [
        r"% Generated by scripts/build_collinearity_latex.py",
        rf"% Source: {os.path.basename(csv_path)}",
        r"% Baseline: pointwise formulation",
        r"% Requires: \usepackage{booktabs}",
        "",
    ]
    for metric in metrics:
        for axis in ("eq_name", "library"):
            out.append(build_table(data, metric, axis, drop_censored))
            out.append("")
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(
        description="Generate LaTeX collinearity penalty tables from a "
                    "diagnostics CSV.",
    )
    parser.add_argument(
        "csv",
        nargs="?",
        default=os.path.join(
            CODE_DIR, "results", "collinearity_diagnostics_spectral.csv",
        ),
        help="Collinearity diagnostics CSV (default: the spectral results).",
    )
    parser.add_argument(
        "--metric", "-m", default="both",
        choices=["max_vif", "mean_vif", "cond_normalised", "both"],
        help="Metric to summarise; 'both' emits max_vif and cond_normalised.",
    )
    parser.add_argument(
        "--drop-censored", action="store_true",
        help="Exclude cells where max_vif sits on the 1e10 solver cap.",
    )
    parser.add_argument("--out", "-o", help="Write to file instead of stdout.")
    args = parser.parse_args()

    metrics = (["max_vif", "cond_normalised"] if args.metric == "both"
               else [args.metric])
    latex = build_document(args.csv, metrics, args.drop_censored)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write(latex + "\n")
        print(f"Written {args.out}")
    else:
        print(latex)


if __name__ == "__main__":
    main()
