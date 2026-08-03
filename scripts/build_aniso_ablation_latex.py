#!/usr/bin/env python3
"""
Generate LaTeX ablation tables for anisotropic noise sweeps.

Outputs copy-pasteable LaTeX ``tabular`` environments for each equation,
noise level, and noise profile, wrapped in ``table`` floats with captions
and labels.

Usage
-----
    python scripts/build_aniso_ablation_latex.py results/sweep_*_aniso*.csv

    # Restrict to one equation
    python scripts/build_aniso_ablation_latex.py results/sweep.csv --equation ac

    # Write to file instead of stdout
    python scripts/build_aniso_ablation_latex.py results/sweep.csv -o tables.tex
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
from equivariant_pde_find.summary_metadata import (  # noqa: E402
    EQUATION_META, get_meta, PROFILE_LABELS,
)
from equivariant_pde_find.markdown_tables import (  # noqa: E402
    build_aniso_ablation_latex,
)

# ── library name normalization (same as build_ablation_latex.py) ──────

_LIB_NAME_NORMALIZE = {
    "multigraph": "multigraphs",
    "mod_multigraph": "multigraphs",
    "pdefind": "standard",
    "polyjet": "polynomial_jet",
    "wsindy": "wsindy",
}

_NOISE_CAPTION = {"Clean": "0% noise", "2%": "2% noise", "5%": "5%"}


def build_latex_document(csv_paths, equation=None):
    """Return a LaTeX document string with anisotropic ablation tables."""
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

    _EQ_ORDER = ["ac", "sh", "pme", "ch"]
    eq_names = [e for e in _EQ_ORDER
                if e in set(r["equation"] for r in all_rows)]
    if equation:
        matches = [e for e in eq_names if e.startswith(equation)]
        if not matches:
            print(f"Equation {equation!r} not found. Available: {eq_names}")
            sys.exit(1)
        eq_names = [matches[0]]

    tables = build_aniso_ablation_latex(all_rows)

    parts = []
    for key, tex in sorted(tables.items()):
        # key format: "{eq}_{profile_type}_{profile_direction}_{noise_label}"
        # e.g. "ac_gaussian_decay_radial_2%"
        parts_key = key.rsplit("_", 3)
        if len(parts_key) < 4:
            continue
        eq_name = parts_key[0]
        pt = parts_key[1]
        pd = parts_key[2]
        noise_label = parts_key[3]

        if eq_name not in eq_names:
            continue
        if not tex:
            continue

        meta = get_meta(eq_name) or {}
        title = meta.get("title", eq_name)
        profile_name = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")

        eq_slug = eq_name.replace("_", "-")
        noise_slug = noise_label.lower().replace("%", "pct")
        profile_slug = f"{pt}_{pd}"
        caption = (
            f"{title} — {profile_name}, "
            f"{_NOISE_CAPTION.get(noise_label, noise_label)}"
        )
        caption = caption.replace("&", r"\&").replace("%", r"\%")
        label = f"tab:{eq_slug}-aniso-{profile_slug}-{noise_slug}"

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
        description="Generate LaTeX anisotropic ablation tables.",
    )
    parser.add_argument(
        "csv", nargs="+",
        help="One or more sweep CSV files to read.",
    )
    parser.add_argument(
        "--equation", "-e",
        help="Restrict to one equation (prefix match).",
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
