"""
Collinearity diagnostics sweep.

Computes max pairwise correlation, max/mean VIF, and condition numbers
(raw and column-normalised) for all library types, formulations, and
derivative backends across benchmark PDEs.

Loads pre-generated .npz data and outputs a CSV of results plus a
structured Markdown summary.

Usage
-----
    # Quick test: Swift-Hohenberg, spectral only, pointwise only
    python scripts/diagnostics_collinearity.py \\
        --equations sh --backends spectral --formulations pointwise

    # Full sweep (all equations, both backends, all formulations)
    OMP_NUM_THREADS=1 python scripts/diagnostics_collinearity.py \\
        --max-snaps 30 \\
        --output results/collinearity_diagnostics.csv \\
        --summary results/collinearity_diagnostics.md
"""

import argparse
import csv
import os
import sys
import time

import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        desc = kwargs.get("desc", "")
        items = list(iterable)
        print(f"{desc} ({len(items)} items)...")
        return items

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives                        # noqa: E402
from equivariant_pde_find.diagnostics import collinearity_report    # noqa: E402
from equivariant_pde_find.summary_metadata import normalize_sweep_row  # noqa: E402
from equivariant_pde_find.library import (                          # noqa: E402
    build_multigraph_library,
    build_mod_multigraph_library,
    build_pdefind_library,
    build_polyjet_library,
    build_wsindy_library,
)

# ── constants ───────────────────────────────────────────────────────────────
RESULTS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
os.makedirs(RESULTS_DIR, exist_ok=True)

SPACE_M = 12
SPACE_P_FIXED = 6

ALL_LIBRARIES = ["multigraph", "mod_multigraph", "pdefind", "wsindy", "polyjet"]
ALL_EQUATIONS = ["ac", "sh", "pme", "ch", "ks_cell", "ks_wchaos"]
ALL_NOISE = ["clean", "noise2pct"]
ALL_FORMULATIONS = ["pointwise", "space_weak_direct", "space_weak_ibp"]
ALL_BACKENDS = ["spectral", "polynomial"]

# Which libraries to skip per equation
EXCLUDED = {
    "pme": ["pdefind"],
    "ch": ["pdefind"],
    "ks_cell": ["pdefind", "wsindy"],
    "ks_wchaos": ["pdefind", "wsindy"],
}

# Equation sort order (roughly: complexity)
EQ_ORDER = {eq: i for i, eq in enumerate(["ac", "ch", "ks_cell", "ks_wchaos", "pme", "sh"])}
LIB_ORDER = {lib: i for i, lib in enumerate(ALL_LIBRARIES)}


# ── data loading ────────────────────────────────────────────────────────────

def load_dataset(npz_path):
    """Load a .npz data file."""
    data = np.load(npz_path)
    snapshots = data["snapshots"]          # (T, N, N)
    times = data["times"]
    Lx = float(data["Lx"])
    Ly = float(data["Ly"])
    N = int(data["N"])
    dt_snap = float(data["dt_snap"])
    obs_noise = float(data["obs_noise"])
    eq_name = str(data.get("eq_name", ""))
    return snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name


# ── library building ────────────────────────────────────────────────────────

def build_library_snapshot(f, Lx, Ly, library, formulation, backend,
                           space_ibp=False, m=SPACE_M, p=SPACE_P_FIXED):
    """Build Theta for a single snapshot.

    Returns (Theta_snap, descriptions) where Theta_snap is (n_points, K).
    """
    weak_kwargs = {}
    if formulation in ("space_weak_direct", "space_weak_ibp"):
        weak_kwargs = {"formulation": "space_weak", "m": m, "p": p, "space_ibp": space_ibp}
    else:
        weak_kwargs = {"formulation": "pointwise"}

    builder_map = {
        "multigraph":     lambda: build_multigraph_library(f, Lx, Ly, **weak_kwargs),
        "mod_multigraph": lambda: build_mod_multigraph_library(f, Lx, Ly, **weak_kwargs),
        "pdefind":        lambda: build_pdefind_library(f, Lx, Ly, max_order=4, max_poly=3, **weak_kwargs),
        "wsindy":         lambda: build_wsindy_library(f, Lx, Ly, max_order=4, max_poly=3, **weak_kwargs),
        "polyjet":        lambda: build_polyjet_library(f, Lx, Ly, max_order=4, max_degree=3, **weak_kwargs),
    }

    return builder_map[library]()


# ── CSV / markdown output ───────────────────────────────────────────────────

CSV_COLUMNS = [
    "eq_name", "noise_level", "library", "formulation", "backend",
    "n_cols", "n_rows", "max_abs_corr", "max_vif", "mean_vif",
    "cond_raw", "cond_normalised", "skipped_cols",
]


def fmt_sci(val):
    """Format a float to 3 significant figures in scientific notation."""
    if np.isnan(val):
        return "nan"
    if np.isinf(val):
        return "inf" if val > 0 else "-inf"
    if val == 0:
        return "0"
    return f"{val:.3g}"


def _write_summary_header(f, formulation, backend):
    """Write a markdown section header for a formulation/backend pair."""
    label_map = {
        "pointwise": "Pointwise",
        "space_weak_direct": "Space-weak Direct",
        "space_weak_ibp": "Space-weak IBP",
    }
    backend_map = {"spectral": "Spectral", "polynomial": "Polynomial"}
    f.write(f"## {label_map[formulation]} — {backend_map[backend]}\n\n")


def _write_noise_block(f, noise_label, rows):
    """Write a markdown table for one noise level."""
    f.write(f"### {noise_label}\n\n")
    header = (
        "| Equation | Library | # cols | max\\|cor\\| | max VIF | mean VIF "
        "| κ_raw | κ_norm |\n"
    )
    sep = (
        "|---|---|---:|---:|---:|---:|---:|\n"
    )
    f.write(header)
    f.write(sep)
    for row in rows:
        eq = row["eq_name"]
        lib = row["library"]
        n_cols = row["n_cols"]
        mc = fmt_sci(row["max_abs_corr"])
        mv = fmt_sci(row["max_vif"])
        av = fmt_sci(row["mean_vif"])
        kr = fmt_sci(row["cond_raw"])
        kn = fmt_sci(row["cond_normalised"])
        f.write(f"| {eq} | {lib} | {n_cols} | {mc} | {mv} | {av} | {kr} | {kn} |\n")
    f.write("\n")


def write_markdown_summary(csv_path, md_path):
    """Read the diagnostics CSV and produce a structured Markdown summary."""
    with open(csv_path, "r", newline="") as fh:
        all_rows = [normalize_sweep_row(r) for r in csv.DictReader(fh)]

    # Parse numeric fields
    for row in all_rows:
        for key in ("max_abs_corr", "max_vif", "mean_vif", "cond_raw", "cond_normalised"):
            val = row.get(key, "")
            try:
                row[key] = float(val)
            except (ValueError, TypeError):
                row[key] = float("nan")
        row["n_cols"] = int(row.get("n_cols", 0))

    with open(md_path, "w") as f:
        f.write("# Collinearity Diagnostics\n\n")
        f.write(
            "Computed on the **full design matrix** (all spatial points, "
            "all available snapshots).  "
            "**VIF** uses the no-intercept formula:\n"
            "  $\\mathrm{VIF}_j = 1/(1-R^2_j)$ where $R^2_j$ is the "
            "$R^2$ from an OLS regression of column $j$ on all other "
            "columns (no intercept).  "
            "Constant columns are excluded from VIF statistics.\n"
            "**κ_raw** = $\\sigma_{\\max}/\\sigma_{\\min}$ of the raw "
            "design matrix; **κ_norm** = same after scaling each column "
            "to unit 2-norm.\n\n"
        )
        f.write(
            f"Weak-form parameters: $m = {SPACE_M}$, $p = {SPACE_P_FIXED}$ "
            "(test-function half-support and polynomial degree).  "
            "IBP moves $\\Delta$ and $\\Delta^2$ derivatives onto the "
            "analytically-known test function where supported by the "
            "library type.\n\n"
        )

        for formulation in ALL_FORMULATIONS:
            for backend in ALL_BACKENDS:
                _write_summary_header(f, formulation, backend)
                for noise_level in ALL_NOISE:
                    noise_label = {"clean": "Clean data", "noise2pct": "2% observation noise"}[noise_level]
                    # Filter rows for this (formulation, backend, noise_level)
                    subset = [
                        r for r in all_rows
                        if r.get("formulation") == formulation
                        and r.get("backend") == backend
                        and r.get("noise_level") == noise_level
                    ]
                    if not subset:
                        f.write(f"### {noise_label}\n\n*No data.*\n\n")
                        continue
                    # Sort: by equation order, then library order
                    subset.sort(key=lambda r: (
                        EQ_ORDER.get(r["eq_name"], 99),
                        LIB_ORDER.get(r["library"], 99),
                    ))
                    _write_noise_block(f, noise_label, subset)

    print(f"Markdown summary written to {md_path}")


# ── main sweep ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Collinearity diagnostics for PDE regression libraries",
    )
    parser.add_argument(
        "--equations", type=str, default="ac,sh,pme,ch,ks_cell,ks_wchaos",
        help="Comma-separated equation keys (default: all six)",
    )
    parser.add_argument(
        "--noise", type=str, default="clean,noise2pct",
        help="Comma-separated noise levels (default: clean,noise2pct)",
    )
    parser.add_argument(
        "--formulations", type=str, default="pointwise,space_weak_direct,space_weak_ibp",
        help="Comma-separated formulations (default: all three)",
    )
    parser.add_argument(
        "--backends", type=str, default="spectral,polynomial",
        help="Comma-separated derivative backends (default: both)",
    )
    parser.add_argument(
        "--max-snaps", type=int, default=0,
        help="Max snapshots to use (0 = all).  Use ~30 for polynomial/CH.",
    )
    parser.add_argument(
        "--output", type=str,
        default=os.path.join(RESULTS_DIR, "collinearity_diagnostics.csv"),
        help="Output CSV path",
    )
    parser.add_argument(
        "--summary", type=str,
        default=os.path.join(RESULTS_DIR, "collinearity_diagnostics.md"),
        help="Output Markdown summary path",
    )
    args = parser.parse_args()

    equations = [e.strip() for e in args.equations.split(",")]
    noise_levels = [n.strip() for n in args.noise.split(",")]
    formulations = [f.strip() for f in args.formulations.split(",")]
    backends = [b.strip() for b in args.backends.split(",")]

    # Build list of all (eq, library, noise, formulation, backend) combos
    combos = []
    for eq in equations:
        excluded = EXCLUDED.get(eq, [])
        libs = [lib for lib in ALL_LIBRARIES if lib not in excluded]
        for lib in libs:
            for noise in noise_levels:
                for form in formulations:
                    for be in backends:
                        combos.append((eq, lib, noise, form, be))

    print(f"Total combinations: {len(combos)}")
    print(f"Equations: {equations}")
    print(f"Noise levels: {noise_levels}")
    print(f"Formulations: {formulations}")
    print(f"Backends: {backends}")
    print(f"Max snaps: {args.max_snaps if args.max_snaps > 0 else 'all'}")
    print()

    # Open CSV
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        for eq, lib, noise, form, be in tqdm(combos, desc="Combinations"):
            # Build .npz path
            if noise == "clean":
                npz_name = f"{eq}_clean.npz"
            else:
                npz_name = f"{eq}_{noise}.npz"
            npz_path = os.path.join(DATA_DIR, npz_name)

            if not os.path.exists(npz_path):
                tqdm.write(f"  SKIP {eq}/{noise}: file not found ({npz_path})")
                continue

            tqdm.write(f"  {eq:10s} {lib:14s} {noise:10s} {form:14s} {be:10s}")

            # Load data
            try:
                snapshots, times, Lx, Ly, N, dt_snap, obs_noise, eq_name = load_dataset(npz_path)
            except Exception as e:
                tqdm.write(f"    ERROR loading: {e}")
                continue

            # Limit snapshots if requested
            if args.max_snaps > 0 and len(snapshots) > args.max_snaps:
                idx = np.linspace(0, len(snapshots) - 1, args.max_snaps, dtype=int)
                snapshots = snapshots[idx]

            # Determine IBP flag
            space_ibp = (form == "space_weak_ibp")

            # Build library across all snapshots (no subsampling)
            t_start = time.time()
            all_theta = []
            descriptions = None

            with derivatives.use_backend(be):
                for s, f in enumerate(snapshots):
                    Theta_s, desc = build_library_snapshot(
                        f, Lx, Ly, lib, form, be, space_ibp=space_ibp,
                    )
                    all_theta.append(Theta_s)
                    if descriptions is None:
                        descriptions = desc

            Theta_full = np.vstack(all_theta).astype(float)
            t_build = time.time() - t_start

            # Compute diagnostics
            report = collinearity_report(Theta_full, descriptions)

            # Write CSV row
            writer.writerow({
                "eq_name": eq,
                "noise_level": noise,
                "library": lib,
                "formulation": form,
                "backend": be,
                "n_cols": report["n_cols"],
                "n_rows": report["n_rows"],
                "max_abs_corr": report["max_abs_corr"],
                "max_vif": report["max_vif"],
                "mean_vif": report["mean_vif"],
                "cond_raw": report["cond_raw"],
                "cond_normalised": report["cond_normalised"],
                "skipped_cols": ",".join(report["skipped_columns"]),
            })

            tqdm.write(
                f"    n_rows={report['n_rows']}  n_cols={report['n_cols']}  "
                f"max|cor|={fmt_sci(report['max_abs_corr'])}  "
                f"maxVIF={fmt_sci(report['max_vif'])}  "
                f"meanVIF={fmt_sci(report['mean_vif'])}  "
                f"κ_raw={fmt_sci(report['cond_raw'])}  "
                f"κ_norm={fmt_sci(report['cond_normalised'])}  "
                f"({t_build:.1f}s)"
            )

    print(f"\nCSV written to {args.output}")

    # Write Markdown summary
    write_markdown_summary(args.output, args.summary)
    print("Done.")


if __name__ == "__main__":
    main()
