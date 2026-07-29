#!/usr/bin/env python3
"""
Generate LaTeX ablation tables from a sweep CSV.

Outputs copy-pasteable LaTeX ``tabular`` environments for each equation
and noise level, wrapped in ``table`` floats with captions and labels.

Usage
-----
    python scripts/build_ablation_latex.py results/sweep_<ts>_<name>.csv

    # Restrict to one equation
    python scripts/build_ablation_latex.py results/sweep.csv --equation ac

    # Write to file instead of stdout
    python scripts/build_ablation_latex.py results/sweep.csv -o tables.tex
"""

import argparse
import csv
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
sys.path.insert(0, SRC_DIR)

from equivariant_pde_find.equations import normalize_eq_name  # noqa: E402
from equivariant_pde_find.summary_metadata import EQUATION_META, get_meta  # noqa: E402
from equivariant_pde_find.markdown_tables import build_ablation_latex  # noqa: E402

# ── library name normalization (same as build_summary.py) ──────────────

_LIB_NAME_NORMALIZE = {
    "multigraph": "multigraphs",
    "mod_multigraph": "multigraphs",
    "pdefind": "standard",
    "polyjet": "polynomial_jet",
    "wsindy": "wsindy",
}


def build_latex_document(csv_paths, equation=None):
    """Return a LaTeX document string with ablation tables for all equations."""
    all_rows = []
    for csv_path in csv_paths:
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                row["equation"] = normalize_eq_name(row["equation"])
                if row.get("library") == "equivariant":
                    row["library"] = "multigraphs"
                row["library_original"] = row["library"]
                row["library"] = _LIB_NAME_NORMALIZE.get(
                    row["library"], row["library"],
                )
                if row.get("algorithm") == "strridge":
                    row["algorithm"] = "stridge"
                if "algorithm" not in row:
                    row["algorithm"] = "stridge"
                if "E" in row and "E_coeff" not in row:
                    row["E_coeff"] = row["E"]
                all_rows.append(row)

    eq_names = sorted(set(r["equation"] for r in all_rows))
    if equation:
        matches = [e for e in eq_names if e.startswith(equation)]
        if not matches:
            print(f"Equation {equation!r} not found. Available: {eq_names}")
            sys.exit(1)
        eq_names = [matches[0]]

    NOISE_CAPTION = {"Clean": "0% noise", "2%": "2% noise", "5%": "5% noise"}

    parts = []
    for eq_name in eq_names:
        meta = get_meta(eq_name) or {}
        title = meta.get("title", eq_name)
        rows = [r for r in all_rows if r["equation"] == eq_name]
        tables = build_ablation_latex(rows)

        for noise_label in ["Clean", "2%", "5%"]:
            tex = tables.get(noise_label, "")
            if not tex:
                continue

            eq_slug = eq_name.replace("_", "-")
            noise_slug = noise_label.lower().replace("%", "pct")
            caption = f"{title} — {NOISE_CAPTION.get(noise_label, noise_label)}"
            # Only need to escape & and % in titles (titles are safe ASCII)
            caption = caption.replace("&", r"\&").replace("%", r"\%")
            label = f"tab:{eq_slug}-{noise_slug}"

            parts.append(
                "\\begin{table}[ht]\n"
                "\\centering\n"
                f"\\caption{{{caption}}}\n"
                f"\\label{{{label}}}\n"
                + tex +
                "\n\\end{table}"
            )

    return "\n\n".join(parts)


def main():
    parser = argparse.ArgumentParser(
        description="Generate LaTeX ablation tables from sweep CSV(s).",
    )
    parser.add_argument(
        "csv", nargs="+",
        help="One or more sweep CSV files to read.",
    )
    parser.add_argument(
        "--equation", "-e",
        help="Restrict to a single equation (prefix match, e.g. 'ac' or 'ks_cell').",
    )
    parser.add_argument(
        "--out", "-o",
        help="Write output to file instead of stdout.",
    )
    args = parser.parse_args()

    latex = build_latex_document(args.csv, args.equation)

    if args.out:
        with open(args.out, "w") as f:
            f.write(latex)
        print(f"Written {args.out}")
    else:
        print(latex)


if __name__ == "__main__":
    main()
