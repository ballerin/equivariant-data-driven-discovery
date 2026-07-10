#!/usr/bin/env python3
"""
Build markdown summaries comparing multigraph vs modified multigraph.

Reads the output of scripts/sweep_multigraph_comparison.py, selects the
best hyperparameters per configuration via lexicographic criterion, and
writes per-equation markdown files.

All tables are built programmatically via ``MarkdownTable`` from
``equivariant_pde_find.markdown_tables``.  Equation metadata and library
labels are centralised in ``equivariant_pde_find.summary_metadata``.

Usage: python scripts/summarize_multigraph_comparison.py results/mg_comparison_*.csv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
SUMMARIES_DIR = os.path.join(CODE_DIR, "summaries")
sys.path.insert(0, SRC_DIR)
os.makedirs(SUMMARIES_DIR, exist_ok=True)

from equivariant_pde_find.summary_metadata import (  # noqa: E402
    EQ_NAMES, EQ_TRUTH, EQ_TERMS, NOISE_LABELS, LIB_NAMES,
)
from equivariant_pde_find.markdown_tables import (  # noqa: E402
    MarkdownTable, best_in_group, compare_bases,
)


# ═══════════════════════════════════════════════════════════════════════════
# Table builders
# ═══════════════════════════════════════════════════════════════════════════

def build_comparison_table(noise_data):
    """Build the best-per-configuration comparison table."""
    t = MarkdownTable(
        ["Noise", "Backend", "m", "IBP", "MG fp/fn", "MG E",
         "Mod fp/fn", "Mod E", "Winner"],
        align=["left"] * 4 + ["center"] * 2 + ["center"] * 2 + ["left"],
    )

    for noise, grp in noise_data:
        ns = NOISE_LABELS.get(noise, noise)
        configs = defaultdict(list)
        for r in grp:
            configs[(r["backend"], r["weak_m"], r["method"])].append(r)

        for (be, mv, mt), _cr in sorted(configs.items()):
            ibp = "yes" if mt == "weak_ibp" else "no"
            mgr = [r for r in configs[(be, mv, mt)]
                   if r["library"] in ("multigraphs",)]
            mdr = [r for r in configs[(be, mv, mt)]
                   if r["library"] == "mod_multigraph"]
            if not mgr or not mdr:
                continue
            bo = best_in_group(mgr)
            bm = best_in_group(mdr)
            w = compare_bases(bo, bm)
            t.add_row(
                ns, be, mv, ibp,
                f"{bo['fp']}/{bo['fn']}", f"{float(bo['E']):.4f}",
                f"{bm['fp']}/{bm['fn']}", f"{float(bm['E']):.4f}",
                w,
            )
    return t


def build_detailed_table(noise_data):
    """Build the detailed per-configuration table."""
    t = MarkdownTable(
        ["Noise", "Backend", "m", "IBP", "Library",
         "fp/fn", "n_act", "E", "Exact?", "d_tol", "cond"],
        align=["left"] * 6 + ["center", "right", "left", "right", "right"],
    )

    for noise, grp in noise_data:
        ns = NOISE_LABELS.get(noise, noise)
        configs = defaultdict(list)
        for r in grp:
            configs[(r["backend"], r["weak_m"], r["method"], r["library"])].append(r)

        for (be, mv, mt, lib), cr in sorted(configs.items()):
            ibp = "yes" if mt == "weak_ibp" else "no"
            best = best_in_group(cr)
            exact = "\u2713" if str(best.get("exact_recovery", "0")) == "1" else ""
            lib_disp = LIB_NAMES.get(lib, lib)
            t.add_row(
                ns, be, mv, ibp, lib_disp,
                f"{best['fp']}/{best['fn']}",
                best["n_active"],
                f"{float(best['E']):.4f}",
                exact,
                best.get("d_tol", ""),
                f"{float(best['cond_theta']):.1f}" if best.get("cond_theta", "nan") != "nan" else "nan",
            )
    return t


def build_truth_table(mg_terms, mod_terms):
    """Build the truth-terms table comparing both bases."""
    t = MarkdownTable(["Library", "True terms"])
    t.add_row("Multigraph", mg_terms)
    t.add_row("Modified MG", mod_terms)
    return t


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="Summarize multigraph comparison sweep.")
    p.add_argument("csv_path")
    p.add_argument("--equation", default=None)
    args = p.parse_args()

    rows = list(csv.DictReader(open(args.csv_path)))
    groups = defaultdict(list)
    for r in rows:
        groups[(r["equation"], r["noise"])].append(r)

    eq_groups = defaultdict(list)
    for (eq, noise), grp in groups.items():
        eq_groups[eq].append((noise, grp))

    for eq_name in sorted(eq_groups):
        if args.equation and args.equation not in eq_name:
            continue
        display = EQ_NAMES.get(eq_name, eq_name)
        truth = EQ_TRUTH.get(eq_name, "")
        mg_terms, mod_terms = EQ_TERMS.get(eq_name, ("?", "?"))

        noise_data = sorted(eq_groups[eq_name], key=lambda x: float(x[0]))

        truth_t = build_truth_table(mg_terms, mod_terms)
        comp_t = build_comparison_table(noise_data)
        detail_t = build_detailed_table(noise_data)

        L = []
        L.append(f"# {display} \u2014 Multigraph vs Modified Multigraph")
        L.append("")
        L.append(f"Equation: {truth}")
        L.append("")
        L.append(truth_t.render())
        L.append("")
        L.append("## Best per configuration (STRidge)")
        L.append("")
        L.append(comp_t.render())
        L.append("")
        L.append("## Detailed results")
        L.append("")
        L.append(detail_t.render())
        L.append("")
        L.append("---")
        L.append(f"*Source: `{os.path.basename(args.csv_path)}`*")
        L.append("")

        out_path = os.path.join(SUMMARIES_DIR, f"mg_comparison_{eq_name}.md")
        with open(out_path, "w") as f:
            f.write("\n".join(L))
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
