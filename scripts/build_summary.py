#!/usr/bin/env python3
"""
Build a per-equation summary markdown file from a sweep CSV file.

Reads the sweep outputs (produced by scripts/sweep.py), selects the
best hyperparameters for each (library, noise, method, backend, algorithm)
group, and writes a markdown file with:

  - Some details about the equation under consideration
  - A representability table (which terms each library can express)
  - A summary table comparing equivariant vs best non-equivariant method
  - Detailed pointwise and weak-form tables

The summary table uses a 5-level colour indicator for whether the
equivariant (multigraph) method wins:

  🟢  equivariant clear win
  🟩  equivariant marginal win
  ⚪  tie
  🟥  non-equivariant marginal win
  🔴  non-equivariant clear win
  (empty)  both methods failed

Usage
-----
    python scripts/build_summary.py results/sweep_<timestamp>_<name>.csv

    # Optional: restrict to one equation
    python scripts/build_summary.py results/sweep.csv --equation PME

The output is written to summaries/<equation>.md.
"""

import argparse
import csv
import os
import sys

# ── paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
SUMMARIES_DIR = os.path.join(CODE_DIR, "summaries")
sys.path.insert(0, SRC_DIR)
os.makedirs(SUMMARIES_DIR, exist_ok=True)

from equivariant_pde_find.equations import normalize_eq_name  # noqa: E402
from equivariant_pde_find.summary_metadata import (  # noqa: E402
    LEGEND, EQUATION_META, get_meta,
)
from equivariant_pde_find.markdown_tables import (  # noqa: E402
    MarkdownTable, _read_log_info,
    build_standard_summary_table, build_pointwise_table, build_weak_table,
    build_oracle_table, build_ablation_table,
)

# Normalize library names from sweep.py output to the legacy names expected
# by markdown_tables.py (which was written for older sweep output formats).
_LIB_NAME_NORMALIZE = {
    "multigraph": "multigraphs",
    "mod_multigraph": "multigraphs",
    "pdefind": "standard",
    "polyjet": "polynomial_jet",
    "wsindy": "wsindy",        # unchanged but explicit
}


# ═══════════════════════════════════════════════════════════════════════════
# Markdown assembly
# ═══════════════════════════════════════════════════════════════════════════

def build_markdown(csv_paths, equation=None):
    """Read one or more sweep CSVs and return the complete markdown string."""
    all_rows = []
    for csv_path in csv_paths:
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                # Normalize legacy column values
                row["equation"] = normalize_eq_name(row["equation"])
                if row.get("library") == "equivariant":
                    row["library"] = "multigraphs"
                # Preserve original library name for display before normalizing
                row["library_original"] = row["library"]
                # Normalize current sweep.py library names to legacy names
                row["library"] = _LIB_NAME_NORMALIZE.get(
                    row["library"], row["library"],
                )
                # Normalize legacy algorithm spelling (strridge → stridge)
                if row.get("algorithm") == "strridge":
                    row["algorithm"] = "stridge"
                # Some sweeps omit algorithm (all rows use STRidge)
                if "algorithm" not in row:
                    row["algorithm"] = "stridge"
                # Legacy CSVs use "E", new CSVs use "E_coeff"
                if "E" in row and "E_coeff" not in row:
                    row["E_coeff"] = row["E"]
                all_rows.append(row)

    # Filter to one equation if requested
    eq_names = sorted(set(r["equation"] for r in all_rows))
    if equation:
        matches = [e for e in eq_names if e.startswith(equation)]
        if not matches:
            print(f"Equation {equation!r} not found. Available: {eq_names}")
            sys.exit(1)
        equation = matches[0]
        rows = [r for r in all_rows if r["equation"] == equation]
    else:
        equation = sorted(eq_names)[0]
        rows = [r for r in all_rows if r["equation"] == equation]

    meta = get_meta(equation)
    if meta is None:
        print(f"No metadata for equation {equation!r}. "
              f"Available: {list(EQUATION_META)}")
        sys.exit(1)

    # ── Representability table (programmatic) ───────────────────────────
    rep_t = MarkdownTable(
        ["Operator", "True coeff", "PDE-FIND", "WSINDy", "Poly jet", "Multigraph"],
    )
    for term, pf, ws, pj, mg in meta["standard_terms"]:
        yn = lambda b: "\U0001f7e2" if b else "\U0001f534"
        rep_t.add_row(f"${term}$", "", yn(pf), yn(ws), yn(pj), yn(mg))

    # ── Equivariant terms table ────────────────────────────────────────
    eq_terms_t = MarkdownTable(["Operator", "True coefficient"])
    for term, coeff in meta["equivariant_terms"]:
        eq_terms_t.add_row(f"${term}$", f"${coeff}$")

    # ── Summary table ─────────────────────────────────────────────────
    summary_t = build_standard_summary_table(rows)

    # ── Ablation tables ───────────────────────────────────────────────
    ablation_tables = build_ablation_table(rows)
    ablation_sections = []
    for noise_label in ["Clean", "2%", "5%"]:
        at = ablation_tables.get(noise_label)
        if at is not None:
            ablation_sections.append(
                f"### {noise_label} (0% noise)\n"
                if noise_label == "Clean"
                else f"### {noise_label} noise\n"
            )
            ablation_sections.append(f"{at.render()}\n")
    ablation_block = (
        "\n## Ablation\n\n"
        "Each row shows the best configuration for a (library, algorithm) "
        "pair, ranked from best to worst (exact recovery, then fn, then fp, "
        "then E_coeff).  "
        "vs. best = (E \u2212 E_best) / E_best \u00d7 100% "
        "(\u2717 if method does not achieve exact recovery).\n\n"
        + "\n".join(ablation_sections)
    ) if ablation_sections else ""

    # ── Pointwise tables ──────────────────────────────────────────────
    pw_spec = build_pointwise_table(rows, "spectral", "Spectral derivatives")
    pw_poly = build_pointwise_table(rows, "polynomial", "Polynomial stencils")

    # ── Weak-form tables ──────────────────────────────────────────────
    wk_spec = build_weak_table(rows, "spectral", "Spectral derivatives")
    wk_poly = build_weak_table(rows, "polynomial", "Polynomial stencils")

    # ── Sources table ──────────────────────────────────────────────────
    src_t = MarkdownTable(
        ["Source", "Wall time", "CPU", "Threads"],
        align=["left", "right", "left", "center"],
    )
    for p in csv_paths:
        info = _read_log_info(p)
        src_t.add_row(
            f"`{os.path.basename(p)}`",
            f"{info.get('total_wall_time_min', '?')} min",
            info.get("cpu_model", "?"),
            info.get("omp_num_threads", "?"),
        )

    # ── Oracle table ───────────────────────────────────────────────────
    oracle_t = build_oracle_table(rows)
    oracle_section = (
        f"\n## Oracle (ground-truth residual)\n\n"
        f"The oracle plugs the true coefficients into "
        f"``evaluate_unified`` (same pipeline as all regression results).  "
        f"Oracle E = 0 always (w ≡ true).  "
        f"E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.\n\n"
        f"{oracle_t.render()}\n"
    ) if oracle_t else ""

    md = f"""# {meta["title"]}

The equation under consideration is

$${meta["latex"]}$$

which can be written in terms of $O(2)$-equivariant terms

{eq_terms_t.render()}

or standard terms (PolyJet basis)

{rep_t.render()}

{meta["notes"]}

{LEGEND}

{summary_t.render()}

{ablation_block}
## Pointwise formulation

### Spectral derivatives

{pw_spec.render()}

### Polynomial stencils

{pw_poly.render()}

## Weak formulation

### Spectral derivatives

{wk_spec.render()}

### Polynomial stencils

{wk_poly.render()}

{oracle_section}---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation {equation}`.

{src_t.render()}
"""
    return md


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Build per-equation summary markdown from a sweep CSV."
    )
    parser.add_argument(
        "csv", nargs="+",
        help="One or more sweep CSV files (e.g. results/sweep_*.csv)"
    )
    parser.add_argument(
        "--equation", "-e", default=None,
        help="Restrict to one equation (default: first found in CSV)."
    )
    parser.add_argument(
        "--out", "-o", default=None,
        help="Output path (default: summaries/<equation>.md)."
    )
    args = parser.parse_args()

    md = build_markdown(args.csv, args.equation)

    # Determine output path
    if args.out:
        out_path = args.out
    else:
        with open(args.csv[0]) as f:
            reader = csv.DictReader(f)
            eq_names = set(r["equation"] for r in reader)
        eq = args.equation or sorted(eq_names)[0]
        out_path = os.path.join(SUMMARIES_DIR, f"{eq}.md")

    with open(out_path, "w") as f:
        f.write(md)

    print(f"Written {out_path}")


if __name__ == "__main__":
    main()
