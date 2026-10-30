#!/usr/bin/env python3
"""Generate a LaTeX table of multicollinearity diagnostics from the CSV file.

Reads ``results/collinearity_diagnostics_spectral.csv`` (produced by
``scripts/diagnostics_collinearity.py``) and prints a LaTeX table to stdout.

Usage::

    conda run -n PDEDisc python scripts/table_collinearity.py > results/table_collinearity.tex
    conda run -n PDEDisc python scripts/table_collinearity.py --equations ac,sh,pme,ch

The table includes clean-data, pointwise-formulation, spectral-backend rows.
By default all equations in the CSV are included; use ``--equations`` to
restrict to a comma-separated list.
"""

import argparse
import csv
import sys
from pathlib import Path

# --- Configuration -----------------------------------------------------------

CSV_PATH = Path(__file__).resolve().parent.parent / "results" / "collinearity_diagnostics_spectral.csv"

# Library display order (all five, though not every equation has all).
LIBRARY_ORDER = ["pdefind", "wsindy", "multigraph", "mod_multigraph", "polyjet"]

LIBRARY_MACROS = {
    "pdefind":        r"\LibFIND",
    "wsindy":         r"\LibWS",
    "multigraph":     r"\LibM",
    "mod_multigraph": r"\LibMM",
    "polyjet":        r"\LibPj",
}

EQUATION_LABELS = {
    "ac":        "AC",
    "sh":        "SH",
    "pme":       "PME",
    "ch":        "CH",
    "ks_cell":   r"KS (cell)",
    "ks_wchaos": r"KS (w.c.)",
}

# --- Helpers ----------------------------------------------------------------

def fmt_val(x: float) -> str:
    """Format a numeric value for the table.

    * Values ≥ 100: integer with comma-separated thousands (no rounding of
      the integer part beyond rounding to the nearest integer).
    * Values 1–99: at most one decimal place; ``.0`` is stripped.
    """
    r = round(x) if x >= 100 else round(x, 1)
    if r >= 100:
        return f"{int(r):,}"
    # 1–99: one decimal place, strip .0
    s = f"{r:.1f}"
    return s[:-2] if s.endswith(".0") else s


def fmt_corr(x: float) -> str:
    """Format max-abs-correlation (always 3 decimal places)."""
    return f"{x:.3f}"


# --- Main -------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--equations", type=str, default=None,
        help="Comma-separated list of equation keys (e.g. ac,sh,pme,ch).  "
             "Default: all equations in the CSV.",
    )
    args = parser.parse_args()

    if not CSV_PATH.exists():
        print(f"ERROR: CSV not found at {CSV_PATH}", file=sys.stderr)
        sys.exit(1)

    with open(CSV_PATH, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    # Filter: clean data, pointwise formulation, spectral backend
    clean_pt = [
        r for r in rows
        if r["noise_level"] == "clean"
        and r["formulation"] == "pointwise"
        and r["backend"] == "spectral"
    ]

    # Index by (eq_name, library)
    index: dict[tuple[str, str], dict] = {}
    for r in clean_pt:
        index[(r["eq_name"], r["library"])] = r

    # Determine which equations are present, preserving CSV order
    seen: set[str] = set()
    eq_order: list[str] = []
    for r in clean_pt:
        eq = r["eq_name"]
        if eq not in seen:
            seen.add(eq)
            eq_order.append(eq)

    # Apply --equations filter if given
    if args.equations is not None:
        wanted = set(args.equations.split(","))
        eq_order = [e for e in eq_order if e in wanted]

    # Build the LaTeX table
    lines: list[str] = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(
        r"\caption{Multicollinearity metrics across equations and libraries "
        r"under zero-noise conditions.  All differential operators are "
        r"computed via the spectral backend.}"
    )
    lines.append(r"\label{tab:collinearity_libraries}")
    lines.append(r"\begin{tabular}{llrrrr}")
    lines.append(r"\toprule")
    lines.append(
        r"Equation & Library & $\max|\mathrm{corr}|$ & "
        r"max VIF & mean VIF & $\kappa_{\mathrm{norm}}$ \\"
    )
    lines.append(r"\midrule")

    for i_eq, eq in enumerate(eq_order):
        libs_present = [lib for lib in LIBRARY_ORDER if (eq, lib) in index]
        n_libs = len(libs_present)

        for j, lib in enumerate(libs_present):
            r = index[(eq, lib)]
            corr = float(r["max_abs_corr"])
            max_vif = float(r["max_vif"])
            mean_vif = float(r["mean_vif"])
            cond_norm = float(r["cond_normalised"])

            eq_cell = (
                rf"\multirow{{{n_libs}}}{{*}}{{{EQUATION_LABELS.get(eq, eq)}}}"
                if j == 0
                else ""
            )
            lib_cell = LIBRARY_MACROS.get(lib, lib)

            lines.append(
                f"  {eq_cell}"
                f"  & ${lib_cell}$"
                f" & {fmt_corr(corr)}"
                f" & {fmt_val(max_vif)}"
                f" & {fmt_val(mean_vif)}"
                f" & {fmt_val(cond_norm)}"
                r" \\"
            )

        if i_eq < len(eq_order) - 1:
            lines.append(r"\midrule")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
