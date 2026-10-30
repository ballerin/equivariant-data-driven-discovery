#!/usr/bin/env python3
"""Summarise the space-time weak sweep with centering as a compared factor.

The sweep runs one mode (``spacetime_weak_ibp``, polynomial backend, periodic
boundary, M&B test functions) over every library that can represent the
equation, with centering off and on.  Centering is a *compared factor*, not a
hyperparameter: selecting over it would be a second oracle, subtler than the
lambda one and just as wrong.  The existing summarisers
(``summarize_best.py``, ``build_summary.py``, ``build_ablation_latex.py``)
have no notion of the ``centering`` column and would group the two arms
together and take the best, which is exactly that mistake -- hence this
script.

Three selections are reported side by side, per equation x noise x library x
solver, separately for each centering arm:

``delta``
    The truth-free rule of new-plan.md 5.1: among the settings whose residual
    is within a factor (1 + delta) of the best, take the one with the fewest
    active terms; ties by lower residual, then lower lambda.  This is the
    rule a practitioner could actually apply.
``oracle``
    ``evaluation.best_of_sweep``'s lexicographic key (exact recovery, then
    fp, fn, E_coeff).  It reads the ground truth, so it is an upper
    reference, reported for continuity with the earlier tables only.
``mean``
    The fraction of all settings in the arm that recover, i.e. what a setting
    drawn at random would give.

Caveats, which belong in any caption built from this output:

* The residual driving the ``delta`` rule is ``E_data``, measured in sample.
  The protocol's held-out block (new-plan.md 4) is not implemented yet, so
  this is the n >> p proxy of 5.2, not a validation residual.
* One realisation per dataset (ic_seed = 42 throughout), so there are no
  error bars over seeds.  McNemar below is over *settings* within a cell,
  which is a paired test of the centering switch, not of the equations.
* ``E_coeff`` is averaged over recovering rows only; it is meaningless where
  the support is wrong.

Usage
-----
    python scripts/summarize_paper.py results/paper_20260918/sweep_*.csv
    python scripts/summarize_paper.py results/paper_20260918/sweep_*.csv \
        --delta 0.5 --out results/paper_20260918/summary.txt \
        --csv results/paper_20260918/summary.csv
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src")))

from equivariant_pde_find.summary_metadata import (  # noqa: E402
    normalize_sweep_row,
)

EQ_ORDER = ["ac", "sh", "ch", "pme", "ks_cell", "ks_wchaos", "ks_chaos"]
# sweep.py writes the numeric level; the labels are accepted too, so
# that a CSV written by an older or relabelled run still orders right.
NOISE_ORDER = ["0.0", "0.02", "0.05", "clean", "noise2pct", "noise5pct"]
SOLVER_ORDER = ["stridge", "rfe", "mstls"]


# ── helpers ─────────────────────────────────────────────────────────────────


def _float(value, default=float("nan")):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def setting_key(row):
    """The hyperparameter setting, i.e. what centering on/off is paired over."""
    return (row["lam"], row["d_tol"], row["mstls_alpha"])


def wilson(k, n, z=1.96):
    """Wilson score interval for k successes out of n."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def mcnemar_exact(b, c):
    """Two-sided exact McNemar p-value for discordant counts (b, c).

    b and c are the numbers of settings that recover under exactly one of the
    two arms.  Under the null the discordant pairs split 50/50, so the test is
    a two-sided binomial test on min(b, c) out of b + c.
    """
    n = b + c
    if n == 0:
        return float("nan")
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2.0 ** n)
    return min(1.0, 2.0 * tail)


def select_delta(rows, delta):
    """Sparsest model whose residual is within (1 + delta) of the best.

    Ties are broken by lower residual, then by lower lambda -- never by
    anything that reads the ground truth.
    """
    usable = [r for r in rows if not math.isnan(r["E_data"])]
    if not usable:
        return None
    floor = min(r["E_data"] for r in usable)
    if math.isinf(floor):
        return None
    cutoff = floor * (1.0 + delta)
    candidates = [r for r in usable if r["E_data"] <= cutoff] or usable
    return min(candidates, key=lambda r: (r["n_active"], r["E_data"],
                                          _float(r["lam"], float("inf"))))


def select_oracle(rows):
    """``evaluation.best_of_sweep``'s lexicographic key, inlined over dicts."""
    if not rows:
        return None
    return min(rows, key=lambda r: (0 if r["exact_recovery"] else 1,
                                    r["fp"], r["fn"], r["E_coeff"]))


def load(paths):
    """Read the sweep CSVs into fitted rows, dropping the oracle reference."""
    rows = []
    for path in paths:
        with open(path, newline="") as f:
            for raw in csv.DictReader(f):
                raw = normalize_sweep_row(raw)
                if raw["algorithm"] == "oracle":
                    continue
                rows.append(dict(
                    raw,
                    n_active=int(_float(raw["n_active"], 0)),
                    fp=_float(raw["fp"], float("inf")),
                    fn=_float(raw["fn"], float("inf")),
                    exact_recovery=bool(int(_float(raw["exact_recovery"], 0))),
                    E_coeff=_float(raw["E_coeff"], float("inf")),
                    E_data=_float(raw["E_data"]),
                    const_active=(None if raw.get("const_active", "") == ""
                                  else bool(int(raw["const_active"]))),
                    intercept_rel=_float(raw.get("intercept_rel", "")),
                ))
    return rows


# ── report ──────────────────────────────────────────────────────────────────


def _mean(values):
    values = [v for v in values if not math.isnan(v) and not math.isinf(v)]
    return sum(values) / len(values) if values else float("nan")


def _fmt(value, spec="{:.3f}", blank="-"):
    if value is None or (isinstance(value, float)
                         and (math.isnan(value) or math.isinf(value))):
        return blank
    return spec.format(value)


def _order(values, preferred):
    known = [v for v in preferred if v in values]
    return known + sorted(v for v in values if v not in preferred)


def build(rows, delta):
    """One record per (equation, noise, library, solver, centering arm)."""
    cells = defaultdict(list)
    for r in rows:
        cells[(r["equation"], r["noise"], r["library"], r["algorithm"],
               r["centering"])].append(r)

    records = []
    for key, arm in sorted(cells.items(), key=lambda kv: str(kv[0])):
        equation, noise, library, algorithm, centering = key
        d = select_delta(arm, delta)
        o = select_oracle(arm)
        recovered = [r for r in arm if r["exact_recovery"]]
        const = [r["const_active"] for r in arm
                 if r["const_active"] is not None]
        records.append(dict(
            equation=equation, noise=noise, library=library,
            algorithm=algorithm, centering=centering,
            n_settings=len(arm),
            frac_recovering=len(recovered) / len(arm) if arm else float("nan"),
            delta_recovery=(None if d is None else int(d["exact_recovery"])),
            delta_fp=(None if d is None else d["fp"]),
            delta_fn=(None if d is None else d["fn"]),
            delta_n_active=(None if d is None else d["n_active"]),
            delta_const=(None if d is None else d["const_active"]),
            delta_lam=(None if d is None else d["lam"]),
            delta_d_tol=(None if d is None else d["d_tol"]),
            delta_alpha=(None if d is None else d["mstls_alpha"]),
            delta_E_coeff=(None if d is None or not d["exact_recovery"]
                           else d["E_coeff"]),
            oracle_recovery=(None if o is None else int(o["exact_recovery"])),
            oracle_fp=(None if o is None else o["fp"]),
            oracle_fn=(None if o is None else o["fn"]),
            E_coeff_recovered=_mean([r["E_coeff"] for r in recovered]),
            frac_const_active=(_mean([float(c) for c in const])
                               if const else float("nan")),
            mean_intercept_rel=_mean([r["intercept_rel"] for r in arm]),
        ))
    return records


def paired(rows):
    """Per (equation, noise, library, solver): centering on vs off, paired."""
    by_cell = defaultdict(dict)
    for r in rows:
        cell = (r["equation"], r["noise"], r["library"], r["algorithm"])
        by_cell[cell][(setting_key(r), r["centering"])] = r

    out = []
    for cell, arms in by_cell.items():
        settings = {k for (k, _) in arms}
        fixed = broken = both = neither = 0
        for s in settings:
            off = arms.get((s, "off"))
            on = arms.get((s, "on"))
            if off is None or on is None:
                continue          # unpaired setting: not comparable, skipped
            if on["exact_recovery"] and not off["exact_recovery"]:
                fixed += 1
            elif off["exact_recovery"] and not on["exact_recovery"]:
                broken += 1
            elif off["exact_recovery"]:
                both += 1
            else:
                neither += 1
        n_paired = fixed + broken + both + neither
        out.append(dict(
            equation=cell[0], noise=cell[1], library=cell[2],
            algorithm=cell[3], n_paired=n_paired,
            fixed=fixed, broken=broken, both=both, neither=neither,
            p_mcnemar=mcnemar_exact(broken, fixed),
        ))
    return out


def write_report(records, pairs, rows, delta, stream):
    w = stream.write
    w(f"Space-time weak sweep summary   delta={delta}\n")
    w(f"{len(rows)} fitted rows, "
      f"{len({(r['equation'], r['noise']) for r in rows})} datasets\n")
    w("Selection rules: delta = sparsest within (1+delta) of the best "
      "in-sample residual (truth-free);\n"
      "oracle = best_of_sweep's lexicographic key (reads the truth, "
      "reference only).\n\n")

    equations = _order({r["equation"] for r in records}, EQ_ORDER)
    noises = _order({r["noise"] for r in records}, NOISE_ORDER)
    solvers = _order({r["algorithm"] for r in records}, SOLVER_ORDER)

    # ── 1. per-cell selections ───────────────────────────────────────────
    w("=" * 100 + "\n1. Selections per cell (centering off | on)\n"
      + "=" * 100 + "\n")
    hdr = (f"{'library':16s}{'solver':9s}"
           f"{'delta-rec':>10s}{'oracle-rec':>11s}{'frac-set':>9s}"
           f"{'E_coeff':>10s}{'const':>7s}")
    for equation in equations:
        for noise in noises:
            sel = [r for r in records
                   if r["equation"] == equation and r["noise"] == noise]
            if not sel:
                continue
            w(f"\n[{equation}  {noise}]\n")
            w(f"{hdr}    {hdr}\n")
            libraries = sorted({r["library"] for r in sel})
            for library in libraries:
                for solver in solvers:
                    line = ""
                    for centering in ("off", "on"):
                        hit = [r for r in sel
                               if r["library"] == library
                               and r["algorithm"] == solver
                               and r["centering"] == centering]
                        if not hit:
                            line += " " * len(hdr) + "    "
                            continue
                        r = hit[0]
                        line += (
                            f"{library:16s}{solver:9s}"
                            f"{_fmt(r['delta_recovery'], '{:d}'):>10s}"
                            f"{_fmt(r['oracle_recovery'], '{:d}'):>11s}"
                            f"{_fmt(r['frac_recovering'], '{:.2f}'):>9s}"
                            f"{_fmt(r['E_coeff_recovered'], '{:.2e}'):>10s}"
                            f"{_fmt(r['frac_const_active'], '{:.2f}'):>7s}"
                            "    ")
                    w(line.rstrip() + "\n")

    # ── 2. paired centering test ─────────────────────────────────────────
    w("\n" + "=" * 100 + "\n2. Centering on vs off, paired over settings "
      "(McNemar, exact)\n" + "=" * 100 + "\n")
    w(f"{'equation':10s}{'noise':10s}{'library':16s}{'solver':9s}"
      f"{'paired':>7s}{'fixed':>7s}{'broken':>7s}{'both':>6s}"
      f"{'neither':>8s}{'p':>9s}\n")
    for p in sorted(pairs, key=lambda r: (EQ_ORDER.index(r["equation"])
                                          if r["equation"] in EQ_ORDER else 99,
                                          r["noise"], r["library"],
                                          r["algorithm"])):
        if p["fixed"] == 0 and p["broken"] == 0:
            continue              # nothing discordant: centering changed nothing
        w(f"{p['equation']:10s}{p['noise']:10s}{p['library']:16s}"
          f"{p['algorithm']:9s}{p['n_paired']:7d}{p['fixed']:7d}"
          f"{p['broken']:7d}{p['both']:6d}{p['neither']:8d}"
          f"{_fmt(p['p_mcnemar'], '{:.4f}'):>9s}\n")

    total_fixed = sum(p["fixed"] for p in pairs)
    total_broken = sum(p["broken"] for p in pairs)
    w(f"\nPooled over all cells: centering fixed {total_fixed}, "
      f"broke {total_broken}, "
      f"exact McNemar p={_fmt(mcnemar_exact(total_broken, total_fixed), '{:.3g}')}\n")
    w("(Settings are not independent within a cell, so the pooled p-value is\n"
      " indicative only -- it is a paired test of the switch, not of the "
      "equations.)\n")

    # ── 3. totals by selection rule ──────────────────────────────────────
    w("\n" + "=" * 100 + "\n3. Recovery by selection rule, over cells\n"
      + "=" * 100 + "\n")
    w(f"{'solver':9s}{'centering':11s}{'cells':>7s}"
      f"{'delta':>9s}{'95% CI':>16s}{'oracle':>9s}{'random':>9s}"
      f"{'const-active':>14s}\n")
    for solver in solvers:
        for centering in ("off", "on"):
            sel = [r for r in records if r["algorithm"] == solver
                   and r["centering"] == centering]
            if not sel:
                continue
            n = len(sel)
            k_delta = sum(1 for r in sel if r["delta_recovery"])
            k_oracle = sum(1 for r in sel if r["oracle_recovery"])
            lo, hi = wilson(k_delta, n)
            w(f"{solver:9s}{centering:11s}{n:7d}"
              f"{k_delta / n:9.3f}{f'[{lo:.2f}, {hi:.2f}]':>16s}"
              f"{k_oracle / n:9.3f}"
              f"{_mean([r['frac_recovering'] for r in sel]):9.3f}"
              f"{_mean([r['frac_const_active'] for r in sel]):14.3f}\n")

    # ── 4. spurious constant ─────────────────────────────────────────────
    w("\n" + "=" * 100 + "\n4. Spurious constant term "
      "(|w0| / mean|f_t|; centering-tol thresholds this)\n" + "=" * 100 + "\n")
    w(f"{'solver':9s}{'centering':11s}{'frac active':>13s}"
      f"{'mean |w0|/mean|f_t|':>22s}\n")
    for solver in solvers:
        for centering in ("off", "on"):
            sel = [r for r in records if r["algorithm"] == solver
                   and r["centering"] == centering]
            if not sel:
                continue
            w(f"{solver:9s}{centering:11s}"
              f"{_mean([r['frac_const_active'] for r in sel]):13.3f}"
              f"{_mean([r['mean_intercept_rel'] for r in sel]):22.3e}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csvs", nargs="+", help="Sweep CSVs to summarise.")
    ap.add_argument("--delta", type=float, default=0.5,
                    help="Residual slack of the truth-free rule (default 0.5, "
                         "the plateau value measured in new-plan.md 5.2).")
    ap.add_argument("--out", default=None,
                    help="Write the report here as well as to stdout.")
    ap.add_argument("--csv", default=None,
                    help="Write the per-cell records as a tidy CSV.")
    args = ap.parse_args()

    rows = load(args.csvs)
    if not rows:
        raise SystemExit("No fitted rows found.")
    records = build(rows, args.delta)
    pairs = paired(rows)

    write_report(records, pairs, rows, args.delta, sys.stdout)
    if args.out:
        with open(args.out, "w") as f:
            write_report(records, pairs, rows, args.delta, f)
        print(f"\nReport → {args.out}")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        print(f"Per-cell records → {args.csv}")


if __name__ == "__main__":
    main()
