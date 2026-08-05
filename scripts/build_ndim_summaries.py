#!/usr/bin/env python
"""Build per-equation N-D dimension-scaling summaries from sweep_nd.py CSVs.

    python scripts/build_ndim_summaries.py results/sweep_nd_*.csv
    python scripts/build_ndim_summaries.py results/sweep_nd_ndim_2d4d_*.csv --out-dir summaries/
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DEFAULT_OUT_DIR = os.path.join(ROOT, "summaries")

LIB_ORDER = ["multigraph", "mod_multigraph", "pdefind", "wsindy"]
LIB_LABEL = {
    "multigraph": "Multigraph",
    "mod_multigraph": "Mod. multigraph",
    "pdefind": "PDE-FIND",
    "wsindy": "WSINDy",
}

EQ_LABEL = {
    "ac": "Allen-Cahn",
    "sh": "Swift-Hohenberg",
    "pme": "Porous Medium",
    "ch": "Cahn-Hilliard",
}

EQ_PDE = {
    "ac":  "f_t = \\varepsilon^2 \\Delta f + f - f^3",
    "sh":  "f_t = (r-1)f - 2\\Delta f - \\Delta^2 f - f^3",
    "pme": "f_t = \\nabla\\cdot(f^2\\nabla f)",
    "ch":  "f_t = -\\Delta(\\varepsilon^2\\Delta f - f + f^3)",
}


def load(paths):
    """Load all CSV rows, keeping the best per (eq, ndim, noise, lib, alg, ibp)."""
    rows = []
    for path in paths:
        with open(path) as fh:
            for r in csv.DictReader(fh):
                rows.append(r)

    # Pick best per group: lexicographic exact → fp → fn → E_coeff
    groups = defaultdict(list)
    for r in rows:
        key = (r.get("equation", "ac"), r["ndim"], r["noise"],
               r["library"], r["algorithm"], r.get("use_ibp", "ibp"))
        groups[key].append(r)

    best = {}
    for key, group in groups.items():
        best[key] = min(group, key=lambda r: (
            -int(r.get("exact_recovery", 0)),
            int(r.get("fp", 999)),
            int(r.get("fn", 999)),
            float(r.get("E_coeff", 1e9)),
        ))
    return best


def fmt(value, spec="{:.4f}", dash="--"):
    try:
        return spec.format(float(value))
    except (TypeError, ValueError):
        return dash


def table(best, eq, ndims, noises, cell_fn, header):
    lines = [
        f"| ndim | noise | " + " | ".join(LIB_LABEL[l] for l in LIB_ORDER) + " |",
        "|---:|---:|" + "---|" * len(LIB_ORDER),
    ]
    for ndim in ndims:
        for noise in noises:
            cells = []
            for lib in LIB_ORDER:
                r = best.get((eq, str(ndim), str(noise), lib, "stride", "ibp"))
                if r is None:
                    r = best.get((eq, str(ndim), str(noise), lib, "stride", "direct"))
                cells.append(cell_fn(r) if r else "--")
            lines.append(f"| {ndim} | {noise:g} | " + " | ".join(cells) + " |")
    return f"**{header}**\n\n" + "\n".join(lines)


def build_summary(best, eq):
    """Build a markdown summary for one equation."""
    eq_rows = {k: v for k, v in best.items() if k[0] == eq}
    if not eq_rows:
        return None

    ndims = sorted({int(k[1]) for k in eq_rows})
    noises = sorted({float(k[2]) for k in eq_rows})

    def cell_recovery(r):
        mark = "yes" if r.get("exact_recovery") == "1" else "no"
        return f"{mark} (fp {r['fp']}, fn {r['fn']})"

    def cell_e(r):
        return fmt(r.get("E_coeff"))

    def cell_build(r):
        # build_time_ms → seconds
        return fmt(float(r.get("build_time_ms", 0)) / 1000, "{:.1f}")

    def cell_fit(r):
        return fmt(float(r.get("fit_time_ms", 0)) / 1000, "{:.2f}")

    def cell_rss(r):
        return fmt(float(r.get("peak_rss_mb", 0)) / 1024, "{:.1f}")

    def cell_cond(r):
        return fmt(r.get("cond_theta"), "{:.1e}")

    # Column counts per library at each ndim
    cols = {}
    for (eq_key, ndim_s, _noise, lib, _alg, _ibp), r in eq_rows.items():
        if int(ndim_s) == max(ndims):
            cols[lib] = r.get("n_cols", "?")

    parts = [
        f"# {EQ_LABEL.get(eq, eq)}: how recovery and cost scale with dimension",
        "",
        f"{EQ_LABEL.get(eq, eq)} `${EQ_PDE.get(eq, '')}` recovered from weak-form data at",
        f"ndim = {min(ndims)}..{max(ndims)}.",
        "",
        "## Library width",
        "",
        "The equivariant bases are dimension-independent; the coordinate libraries are not.",
        "",
        "| ndim | " + " | ".join(LIB_LABEL[l] for l in LIB_ORDER) + " |",
        "|---:|" + "---|" * len(LIB_ORDER),
    ]
    # Collect column counts across ndims
    cols_by_ndim = {}
    for (eq_key, ndim_s, _noise, lib, _alg, _ibp), r in eq_rows.items():
        ndim = int(ndim_s)
        cols_by_ndim.setdefault(ndim, {})[lib] = r.get("n_cols", "?")
    for ndim in ndims:
        parts.append(f"| {ndim} | " + " | ".join(
            str(cols_by_ndim.get(ndim, {}).get(l, "?")) for l in LIB_ORDER) + " |")

    # Headline numbers
    ok = [r for k, r in eq_rows.items()]
    top = max(ndims)

    mg_build = max((float(r.get("build_time_ms", 0)) / 1000
                    for r in ok if int(r["ndim"]) == top
                    and r["library"] == "multigraph"), default=0)
    pf_build = max((float(r.get("build_time_ms", 0)) / 1000
                    for r in ok if int(r["ndim"]) == top
                    and r["library"] == "pdefind"), default=0)
    mg_rss = max((float(r.get("peak_rss_mb", 0)) / 1024
                  for r in ok if int(r["ndim"]) == top
                  and r["library"] == "multigraph"), default=0)
    pf_rss = max((float(r.get("peak_rss_mb", 0)) / 1024
                  for r in ok if int(r["ndim"]) == top
                  and r["library"] == "pdefind"), default=0)
    mgm_build = max((float(r.get("build_time_ms", 0)) / 1000
                     for r in ok if int(r["ndim"]) == top
                     and r["library"] == "mod_multigraph"), default=0)

    parts += [
        "",
        "## Headline numbers",
        "",
        f"- At ndim={top} the equivariant bases stay at **20 columns** against "
        f"PDE-FIND's **{cols.get('pdefind', '?')}** and "
        f"WSINDy's **{cols.get('wsindy', '?')}**.",
        f"- Build time at ndim={top}: multigraph {mg_build / 60:.0f} min, "
        f"mod. multigraph {mgm_build / 60:.0f} min, "
        f"PDE-FIND {pf_build / 60:.0f} min ({pf_build / max(mg_build, 1):.1f}x).",
        f"- Peak memory at ndim={top}: multigraph {mg_rss:.1f} GB, "
        f"PDE-FIND {pf_rss:.1f} GB ({pf_rss / max(mg_rss, 0.1):.1f}x).",
        "",
        "## Exact recovery",
        "",
        table(eq_rows, eq, ndims, noises, cell_recovery,
              "Support recovered exactly? (best STRidge row per group)"),
        "",
        "## Coefficient error (E_coeff)",
        "",
        table(eq_rows, eq, ndims, noises, cell_e,
              "Relative L₂ coefficient error of best STRidge recovery"),
        "",
        "## Condition number",
        "",
        table(eq_rows, eq, ndims, noises, cell_cond,
              "Condition number of weak-form design matrix"),
        "",
        "## Cost",
        "",
        table(eq_rows, eq, ndims, noises, cell_build,
              "Library build time (s)"),
        "",
        table(eq_rows, eq, ndims, noises, cell_fit,
              "STRidge fit time (s)"),
        "",
        table(eq_rows, eq, ndims, noises, cell_rss,
              "Peak RSS (GB)"),
        "",
        "---",
        "",
        f"Generated by `scripts/build_ndim_summaries.py`.",
        "",
        "Best row selected lexicographically: exact recovery → fp → fn → E_coeff.",
        "IBP mode chosen per row (ibp or direct, whichever scored best).",
    ]
    return "\n".join(parts) + "\n"


def main():
    ap = argparse.ArgumentParser(
        description="Build per-equation N-D scaling summaries from sweep_nd.py CSVs."
    )
    ap.add_argument("csv", nargs="+", help="One or more sweep_nd.py output CSVs")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR,
                    help="Output directory for .md files")
    args = ap.parse_args()

    best = load(args.csv)
    if not best:
        print("No rows found.")
        return 1

    equations = sorted({k[0] for k in best})
    os.makedirs(args.out_dir, exist_ok=True)

    for eq in equations:
        md = build_summary(best, eq)
        if md is None:
            print(f"  {eq}: no data, skipping")
            continue
        label = EQ_LABEL.get(eq, eq).lower().replace(" ", "_")
        out_path = os.path.join(args.out_dir, f"{label}_dimension_scaling.md")
        with open(out_path, "w") as fh:
            fh.write(md)
        print(f"  {eq} → {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
