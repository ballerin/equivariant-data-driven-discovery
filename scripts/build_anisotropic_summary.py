#!/usr/bin/env python3
"""
Build a per-equation anisotropic-noise summary markdown from a sweep CSV.

Reads a sweep CSV containing anisotropic (space-dependent) noise datasets,
selects the best hyperparameters for each (profile, library, noise, method,
backend, algorithm) group, and writes a markdown file comparing equivariant vs
standard (and other) libraries under anisotropic noise.

All tables are built programmatically via ``MarkdownTable`` from
``equivariant_pde_find.markdown_tables``.  Equation metadata is centralised in
``equivariant_pde_find.summary_metadata``.

Usage
-----
    python scripts/build_anisotropic_summary.py results/sweep_<ts>_aniso.csv
    python scripts/build_anisotropic_summary.py results/sweep.csv -e allen_cahn

The output is written to ``summaries/<equation>_anisotropic.md``.
"""

import argparse
import csv
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
SUMMARIES_DIR = os.path.join(CODE_DIR, "summaries")
sys.path.insert(0, SRC_DIR)
os.makedirs(SUMMARIES_DIR, exist_ok=True)

from equivariant_pde_find.summary_metadata import (  # noqa: E402
    LEGEND, LIB_NAMES, METHOD_MAP, PROFILE_LABELS, PROFILE_ORDER,
    EQUATION_META, get_meta,
)
from equivariant_pde_find.markdown_tables import (  # noqa: E402
    MarkdownTable, _read_log_footer,
    build_anisotropic_summary_table, build_anisotropic_detailed_table,
)


# ═══════════════════════════════════════════════════════════════════════════
# Markdown assembly
# ═══════════════════════════════════════════════════════════════════════════

def build_markdown(csv_path, equation=None):
    with open(csv_path) as f:
        all_rows = list(csv.DictReader(f))

    eq_names = sorted(set(r["equation"] for r in all_rows))
    if equation:
        matches = [e for e in eq_names if e.startswith(equation)]
        if not matches:
            print(f"Equation {equation!r} not found. Available: {eq_names}")
            sys.exit(1)
        eq_names = matches

    methods = sorted(set(r["method"] for r in all_rows))
    algorithms = sorted(set(r["algorithm"] for r in all_rows))
    libraries = sorted(set(r["library"] for r in all_rows))
    backends = sorted(set(r["backend"] for r in all_rows))
    noise_levels = sorted(set(r["noise"] for r in all_rows), key=float)
    profiles_found = sorted(set(
        (r.get("profile_type", ""), r.get("profile_direction", ""))
        for r in all_rows if r.get("profile_type", "")
    ))

    sweep_name = os.path.basename(csv_path)

    def _build_one(eq_name):
        eq_rows = [r for r in all_rows if r["equation"] == eq_name]
        meta = get_meta(eq_name)
        if meta is None:
            print(f"No metadata for {eq_name!r}. Skipping.")
            return None

        # ── Representability table ───────────────────────────────────────
        rep_t = MarkdownTable(
            ["Operator", "True coeff", "PDE-FIND", "WSINDy", "Poly jet", "Multigraph"],
        )
        for term, pf, ws, pj, mg in meta["standard_terms"]:
            yn = lambda b: "\U0001f7e2" if b else "\U0001f534"
            rep_t.add_row(f"${term}$", "", yn(pf), yn(ws), yn(pj), yn(mg))

        # ── Equivariant terms table ──────────────────────────────────────
        eq_terms_t = MarkdownTable(["Operator", "True coefficient"])
        for term, coeff in meta["equivariant_terms"]:
            eq_terms_t.add_row(f"${term}$", f"${coeff}$")

        methods_str = ", ".join(METHOD_MAP.get(m, m) for m in methods)
        algs_str = ", ".join(algorithms)
        libs_str = ", ".join(LIB_NAMES.get(l, l) for l in libraries)
        profiles_str = ", ".join(
            PROFILE_LABELS.get(p, f"{p[0]}/{p[1]}") for p in profiles_found
        )

        n_libs = len(libraries)
        lib_plural = "ies" if n_libs > 1 else "y"

        # ── Summary table ────────────────────────────────────────────────
        summary_t = build_anisotropic_summary_table(eq_rows)

        # ── Detailed tables ──────────────────────────────────────────────
        pw_spec = build_anisotropic_detailed_table(
            eq_rows, "pointwise", "pointwise", "spectral", "Spectral derivatives",
        )
        pw_poly = build_anisotropic_detailed_table(
            eq_rows, "pointwise", "pointwise", "polynomial", "Polynomial stencils",
        )

        md = f"""# {meta["title"]} \u2014 Anisotropic noise

The equation under consideration is

$${meta["latex"]}$$

which can be written in terms of $O(2)$-equivariant terms

{eq_terms_t.render()}

or standard terms (PolyJet basis)

{rep_t.render()}

All libraries are evaluated in the common polynomial-jet basis (105 monomials,
max order 4, max degree 3).  fp/fn/E are computed against the same ground truth
for every library.  {meta["notes"]}

**Noise model**: Space-dependent observation noise
$\\sigma \\cdot P(x,y) \\cdot \\mathcal{{N}}(0,1)$ where
$\\sigma = \\mathrm{{obs\\_noise}} \\cdot \\mathrm{{std}}(f)$ and $P(x,y)$ is a
normalised spatial profile ($\\mathrm{{mean}}(P^2)=1$, preserving total noise
power relative to isotropic).

| Profile | $P(x,y)$ | Symmetry |
|---|---|---|
| Radial decay | $\\propto \\exp(-r^2/2s^2)$, $r=\\sqrt{{x^2+y^2}}$, $s=L/4$ | O(2) |
| Radial growth | $\\propto 1 + (1 - \\exp(-r^2/2s^2))$ | O(2) |
| Axial decay | $\\propto \\exp(-x^2/2s^2)$ | Reflection only |
| Axial growth | $\\propto 1 + (1 - \\exp(-x^2/2s^2))$ | Reflection only |

Sweep configuration: {n_libs} librar{lib_plural} ({libs_str}),
{len(algorithms)} algorithm(s) ({algs_str}),
{len(methods)} method(s) ({methods_str}),
{len(backends)} backend(s) ({', '.join(backends)}),
{len(noise_levels)} noise level(s),
{len(profiles_found)} noise profile(s) ({profiles_str}).

{LEGEND}

{summary_t.render()}

## Pointwise formulation

### Spectral derivatives

{pw_spec.render()}

### Polynomial stencils

{pw_poly.render()}

---
*Source: `{sweep_name}`*
{_read_log_footer(csv_path)}
"""
        return md

    for eq in eq_names:
        md = _build_one(eq)
        if md is None:
            continue
        out_path = os.path.join(SUMMARIES_DIR, f"{eq}_anisotropic.md")
        with open(out_path, "w") as f:
            f.write(md)
        print(f"Written {out_path}")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Build anisotropic-noise summary markdown from a sweep CSV."
    )
    parser.add_argument("csv", help="Path to sweep CSV.")
    parser.add_argument("--equation", "-e", default=None,
                        help="Restrict to one equation (default: all found).")
    parser.add_argument("--out", "-o", default=None,
                        help="Output path (default: summaries/<eq>_anisotropic.md).")
    args = parser.parse_args()

    if args.out:
        # Single-file mode
        md = build_markdown(args.csv, args.equation)
        with open(args.out, "w") as f:
            f.write(md)
        print(f"Written {args.out}")
    else:
        build_markdown(args.csv, args.equation)


if __name__ == "__main__":
    main()
