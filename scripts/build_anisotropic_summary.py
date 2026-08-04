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
    _collect_aniso_ablation_entries, best_in_group, fmt_val,
)
from equivariant_pde_find.summary_metadata import (  # noqa: E402
    ALGORITHM_LABELS, BACKEND_LABELS, NOISE_LABELS,
)

# Map internal library names (from sweep CSV) to display names
_LIB_NAME_NORMALIZE = {
    "multigraph": "multigraphs",
    "mod_multigraph": "multigraphs",
    "pdefind": "standard",
    "polyjet": "polynomial_jet",
    "wsindy": "wsindy",
}


# ═══════════════════════════════════════════════════════════════════════════
# Markdown assembly
# ═══════════════════════════════════════════════════════════════════════════


def _build_aniso_detailed_table(rows, method, backend):
    """Build a detailed per-library table for one (method, backend).

    Dynamically discovers libraries, noise levels, algorithms, and
    profiles present in *rows*.  Includes a ``weak_m`` column when
    weak-form rows are present.
    """
    from equivariant_pde_find.summary_metadata import PROFILE_ORDER, PROFILE_LABELS

    # Discover what is actually in the data
    _lib_order = ["standard", "wsindy", "polynomial_jet", "multigraphs"]
    libs_in_data = sorted(set(
        r["library"] for r in rows
        if r["method"] == method and r["backend"] == backend
    ), key=lambda x: _lib_order.index(x) if x in _lib_order else 99)

    noises_in_data = sorted(set(
        str(r["noise"]) for r in rows
    ), key=lambda x: float(x) if x.replace(".", "").replace("-", "").isdigit() else 999)

    algs_in_data = sorted(set(
        r["algorithm"] for r in rows
        if r["method"] == method and r["backend"] == backend
    ))

    profiles_in_data = sorted(set(
        (r.get("profile_type", ""), r.get("profile_direction", ""))
        for r in rows
        if r.get("profile_type", "") and r["method"] == method
    ))

    # Check if weak_m is relevant
    has_weak_m = any(
        str(r.get("weak_m", "")).strip() != ""
        for r in rows if r["method"] == method
    )

    cols = ["Library", "Noise", "Algorithm"]
    if has_weak_m:
        cols.append("m")
    cols.extend(["d_tol/lam", "Profile", "fp/fn", "n_active",
                  "E_coeff", "Time (s)", "Exact?"])

    t = MarkdownTable(cols)

    for lib in libs_in_data:
        lib_label = LIB_NAMES.get(lib, lib)
        for nl_raw in noises_in_data:
            nl_display = NOISE_LABELS.get(nl_raw, f"{float(nl_raw)*100:.0f}%")
            for alg in algs_in_data:
                alg_label = ALGORITHM_LABELS.get(alg, alg)
                for pt, pd in (profiles_in_data if profiles_in_data else PROFILE_ORDER):
                    plabel = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")
                    g = [
                        r for r in rows
                        if r.get("profile_type", "") == pt
                        and r.get("profile_direction", "") == pd
                        and r["library"] == lib
                        and str(r["noise"]) == str(nl_raw)
                        and r["method"] == method
                        and r["backend"] == backend
                        and r["algorithm"] == alg
                        and r.get("fp", "nan") != "nan"
                    ]
                    if not g:
                        continue
                    b = best_in_group(g)
                    fpf = f"{b['fp']}/{b['fn']}"
                    ex = "\U0001f7e2" if str(b.get("exact_recovery", "0")) == "1" else "\U0001f534"
                    E = fmt_val(b["E_coeff"])
                    tim = fmt_val(float(b.get("time_ms", 0)) / 1000, 1)
                    wm = str(b.get("weak_m", "")) if has_weak_m else None

                    if alg == "rfe":
                        param = str(b.get("lam", ""))
                    else:
                        param = str(b.get("d_tol", ""))

                    row_vals = [lib_label, nl_display, alg_label]
                    if has_weak_m:
                        row_vals.append(wm)
                    row_vals.extend([param, plabel, fpf, str(b["n_active"]), E, tim, ex])
                    t.add_row(*row_vals)

    return t if len(t._rows) > 0 else None


def _render_aniso_ablation_md(ablation_entries):
    """Render anisotropic ablation entries as markdown sections.

    *ablation_entries* is the dict returned by
    ``_collect_aniso_ablation_entries``, keyed by
    ``(noise_label, profile_type, profile_direction)`` tuples.
    Each value is a list of entry dicts with keys:
    label, conf, fpf, exact, E_str, time_str, vs_str.
    """
    from equivariant_pde_find.summary_metadata import PROFILE_LABELS

    parts = []

    # Group by noise label first, then by profile
    # Keys are (noise_label, pt, pd) tuples
    by_noise = {}
    for key, entries in ablation_entries.items():
        if isinstance(key, tuple) and len(key) == 3:
            noise_label, pt, pd = key
        else:
            continue
        if noise_label not in by_noise:
            by_noise[noise_label] = []
        by_noise[noise_label].append((pt, pd, entries))

    for noise_label in sorted(by_noise):
        nl_display = NOISE_LABELS.get(noise_label, noise_label)

        for pt, pd, entries in by_noise[noise_label]:
            if not entries:
                continue
            plabel = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")

            t = MarkdownTable([
                "Rank", "Library + Algorithm", "Configuration",
                "fp/fn", "E_coeff", "Time (s)", "vs. best",
            ])

            for i, entry in enumerate(entries):
                rank = i + 1
                label = entry.get("label", "")
                conf = entry.get("conf", "")
                fpf = entry.get("fpf", "")
                exact = entry.get("exact", False)
                E_str = entry.get("E_str", "")
                time_str = entry.get("time_str", "")
                vs_str = entry.get("vs_str", "")

                # Add checkmark for exact recovery
                if exact:
                    fpf += " \u2713"

                t.add_row(
                    str(rank), label, conf,
                    fpf, E_str, time_str, vs_str,
                )

            parts.append(
                f"### {nl_display} \u2014 {plabel}\n\n{t.render()}\n"
            )

    return "\n".join(parts) if parts else "*No ablation data available.*"

def build_markdown(csv_path, equation=None):
    with open(csv_path) as f:
        all_rows = list(csv.DictReader(f))

    # Normalise library names (internal -> display)
    for row in all_rows:
        row["library_original"] = row.get("library", "")
        row["library"] = _LIB_NAME_NORMALIZE.get(
            row.get("library", ""), row.get("library", ""),
        )

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
        algs_str = ", ".join(ALGORITHM_LABELS.get(a, a) for a in algorithms)
        libs_str = ", ".join(LIB_NAMES.get(l, l) for l in libraries)
        profiles_str = ", ".join(
            PROFILE_LABELS.get(p, f"{p[0]}/{p[1]}") for p in profiles_found
        )

        n_libs = len(libraries)
        lib_plural = "ies" if n_libs > 1 else "y"

        # Summary table
        summary_t = build_anisotropic_summary_table(eq_rows)

        # Ablation table
        ablation_entries = _collect_aniso_ablation_entries(eq_rows)
        ablation_tables_md = _render_aniso_ablation_md(ablation_entries)

        # Detailed per-method tables (dynamic: detects what is in the CSV)
        detailed_sections = []
        present_methods = sorted(set(r["method"] for r in eq_rows))
        for method in present_methods:
            method_label = METHOD_MAP.get(method, method)
            for be in sorted(backends):
                be_label = BACKEND_LABELS.get(be, be)
                t = _build_aniso_detailed_table(eq_rows, method, be)
                if t is not None:
                    detailed_sections.append(
                        f"### {method_label} \u2014 {be_label}\n\n"
                        f"{t.render()}\n"
                    )

        detailed_md = "\n".join(detailed_sections) if detailed_sections else (
            "*No detailed tables available.*"
        )

        md = f"""# {meta["title"]} \u2014 Anisotropic noise

The equation under consideration is

$${meta["latex"]}$$

which can be written in terms of $O(2)$-equivariant terms

{eq_terms_t.render()}

or standard terms (PolyJet basis)

{rep_t.render()}

All libraries are evaluated in the common polynomial-jet basis (max order 4,
max degree 3).  fp/fn/E are computed against the same ground truth for every
library.  {meta["notes"]}

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
{len(backends)} backend(s) ({', '.join(BACKEND_LABELS.get(b, b) for b in backends)}),
{len(noise_levels)} noise level(s),
{len(profiles_found)} noise profile(s) ({profiles_str}).

{LEGEND}

{summary_t.render()}

## Ablation

{ablation_tables_md}

## Detailed results

{detailed_md}

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

    # Always use default path (summaries/<eq>_anisotropic.md).
    # The --out flag writes a single file but is not yet implemented
    # for the per-equation loop; use the default path instead.
    build_markdown(args.csv, args.equation)


if __name__ == "__main__":
    main()
