#!/usr/bin/env python3
"""
Build markdown summaries comparing multigraph vs modified multigraph.

Reads the output of scripts/sweep_multigraph_comparison.py, selects the
best hyperparameters per configuration via lexicographic criterion, and
writes per-equation markdown files.

Usage: python scripts/summarize_multigraph_comparison.py results/mg_comparison_*.csv
"""

import argparse, csv, os, sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SUMMARIES_DIR = os.path.join(CODE_DIR, "summaries")
os.makedirs(SUMMARIES_DIR, exist_ok=True)

EQ_NAMES = {
    "allen_cahn": "Allen-Cahn", "swift_hohenberg": "Swift-Hohenberg",
    "cahn_hilliard": "Cahn-Hilliard", "porous_medium": "Porous Medium",
    "ks_cellular": "KS (cellular)", "ks_weak_chaotic": "KS (weak-chaotic)",
    "ks_chaotic": "KS (chaotic)",
}
EQ_TRUTH = {
    "allen_cahn": r"$f_t = \varepsilon^2 \Delta f + f - f^3$",
    "swift_hohenberg": r"$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3$",
    "cahn_hilliard": r"$f_t = -\Delta^2 f + \Delta(f^3 - f)$",
    "porous_medium": r"$f_t = \Delta(f^3)/3$ (mod) = $f^2\Delta f + 2f\|\nabla f\|^2$ (mg)",
    "ks_cellular": r"$f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$",
    "ks_weak_chaotic": r"$f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$",
    "ks_chaotic": r"$f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$",
}
EQ_TERMS = {
    "allen_cahn": ("f, f\u00b3, \u0394f", "f, f\u00b3, \u0394f"),
    "swift_hohenberg": ("f, f\u00b3, \u0394f, \u0394\u00b2f", "f, f\u00b3, \u0394f, \u0394\u00b2f"),
    "ks_cellular": ("\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f", "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f"),
    "ks_weak_chaotic": ("\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f", "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f"),
    "ks_chaotic": ("\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f", "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f"),
    "cahn_hilliard": ("\u0394f, f\u00b2\u00b7\u0394f, f\u00b7\u2016\u2207f\u2016\u00b2, \u0394\u00b2f", "\u0394f, \u0394f\u00b3, \u0394\u00b2f"),
    "porous_medium": ("f\u00b2\u00b7\u0394f, f\u00b7\u2016\u2207f\u2016\u00b2", "\u0394f\u00b3"),
}
NOISE_LABELS = {"0.0": "clean", "0.02": "2%", "0.05": "5%",
                "0.01": "1%", "0.005": "0.5%", "0.001": "0.1%"}
LIB_NAMES = {"multigraph": "Multigraph", "equivariant": "Multigraph",
             "mod_multigraph": "Modified MG"}


def select_best(rows):
    best = rows[0]
    for r in rows[1:]:
        a = (int(r["exact_recovery"]), -int(r["fp"]), -int(r["fn"]), -float(r["E"]))
        b = (int(best["exact_recovery"]), -int(best["fp"]), -int(best["fn"]), -float(best["E"]))
        if a > b:
            best = r
    return best


def compare(o, m):
    eo, em = int(o["exact_recovery"]), int(m["exact_recovery"])
    if em and not eo: return "\U0001f7e2"
    if eo and not em: return "\U0001f534"
    so, sm = int(o["fp"]) + int(o["fn"]), int(m["fp"]) + int(m["fn"])
    if sm < so: return "\U0001f7e2"
    if so < sm: return "\U0001f534"
    if sm == so and float(m["E"]) < float(o["E"]): return "\U0001f7e9"
    if sm == so and float(o["E"]) < float(m["E"]): return "\U0001f7e5"
    return "\u26aa"


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

        L = []
        L.append(f"# {display} \u2014 Multigraph vs Modified Multigraph")
        L.append("")
        L.append(f"Equation: {truth}")
        L.append("")
        L.append("| Library | True terms |")
        L.append("|---|---|")
        L.append(f"| Multigraph | {mg_terms} |")
        L.append(f"| Modified MG | {mod_terms} |")
        L.append("")

        noise_data = sorted(eq_groups[eq_name], key=lambda x: float(x[0]))

        L.append("## Best per configuration (STRidge)")
        L.append("")
        L.append("| Noise | Backend | m | IBP | MG fp/fn | MG E | Mod fp/fn | Mod E | Winner |")
        L.append("|---|---|---|---|---|---|---|---|---|")

        for noise, grp in noise_data:
            ns = NOISE_LABELS.get(noise, noise)
            configs = defaultdict(list)
            for r in grp:
                configs[(r["backend"], r["weak_m"], r["method"])].append(r)
            for (be, mv, mt), _cr in sorted(configs.items()):
                ibp = "yes" if mt == "weak_ibp" else "no"
                mgr = [r for r in configs[(be, mv, mt)] if r["library"] in ("multigraph", "equivariant")]
                mdr = [r for r in configs[(be, mv, mt)] if r["library"] == "mod_multigraph"]
                if not mgr or not mdr: continue
                bo, bm = select_best(mgr), select_best(mdr)
                w = compare(bo, bm)
                L.append(f"| {ns} | {be} | {mv} | {ibp} | {bo['fp']}/{bo['fn']} | {float(bo['E']):.4f} | {bm['fp']}/{bm['fn']} | {float(bm['E']):.4f} | {w} |")
        L.append("")

        L.append("## Detailed results")
        L.append("")
        L.append("| Noise | Backend | m | IBP | Library | fp/fn | n_act | E | Exact? | d_tol | cond |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")

        for noise, grp in noise_data:
            ns = NOISE_LABELS.get(noise, noise)
            configs = defaultdict(list)
            for r in grp:
                configs[(r["backend"], r["weak_m"], r["method"], r["library"])].append(r)
            for (be, mv, mt, lib), cr in sorted(configs.items()):
                ibp = "yes" if mt == "weak_ibp" else "no"
                best = select_best(cr)
                exact = "\u2713" if best["exact_recovery"] == "1" else ""
                lib_disp = LIB_NAMES.get(lib, lib)
                L.append(f"| {ns} | {be} | {mv} | {ibp} | {lib_disp} | {best['fp']}/{best['fn']} | {best['n_active']} | {float(best['E']):.4f} | {exact} | {best['d_tol']} | {float(best['cond_theta']):.1f} |")

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
