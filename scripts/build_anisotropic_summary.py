#!/usr/bin/env python3
"""
Build a per-equation anisotropic-noise summary markdown from a sweep CSV.

Reads a sweep CSV containing anisotropic (space-dependent) noise datasets,
selects the best hyperparameters for each (profile, library, noise, method,
backend, algorithm) group, and writes a markdown file comparing equivariant vs
standard (and other) libraries under anisotropic noise.

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
SUMMARIES_DIR = os.path.join(CODE_DIR, "summaries")
os.makedirs(SUMMARIES_DIR, exist_ok=True)

LIB_MAP = {
    "equivariant": "Multigraphs",
    "standard": "PDE-FIND",
    "polynomial_jet": "PolyJet",
    "wsindy": "WSINDy",
    "wsindy_proj": "WSINDy \u2192 O(2)",
    "mwsindy": "MWSINDy",
}

METHOD_MAP = {
    "pointwise": "pointwise",
    "weak_direct": "weak (no IBP)",
    "weak_ibp": "weak (IBP)",
}

PROFILE_LABELS = {
    ("gaussian_decay", "radial"): "Radial decay",
    ("gaussian_growth", "radial"): "Radial growth",
    ("gaussian_decay", "axial"): "Axial decay",
    ("gaussian_growth", "axial"): "Axial growth",
}

PROFILE_ORDER = [
    ("gaussian_decay", "radial"),
    ("gaussian_growth", "radial"),
    ("gaussian_decay", "axial"),
    ("gaussian_growth", "axial"),
]

EQUATION_META = {
    "allen_cahn": {
        "title": "Allen-Cahn equation",
        "latex": (
            "f_t = f - f^3 + \\varepsilon^2\\Delta f, "
            "\\qquad \\varepsilon = 0.15"
        ),
        "equivariant_terms": [
            ("f", "+1.0"),
            ("f^3", "-1.0"),
            ("\\Delta f", "+\\varepsilon^2 = +0.0225"),
        ],
        "standard_terms": [
            ("f", True, True, True, True),
            ("f^3", True, True, True, True),
            ("\\partial_{xx} f", True, True, True, True),
            ("\\partial_{yy} f", True, True, True, True),
        ],
        "notes": (
            "Allen-Cahn contains only pure monomials ($f$, $f^3$) "
            "and pure second derivatives ($f_{xx}$, $f_{yy}$).  "
            "All four libraries can represent every term, so the "
            "structural advantage of the equivariant basis is "
            "less pronounced here; differences come from "
            "numerical conditioning, not representability."
        ),
    },
    "swift_hohenberg": {
        "title": "Swift-Hohenberg equation",
        "latex": (
            "f_t = (r-1)f - 2\\Delta f - \\Delta^2 f - f^3, "
            "\\qquad r = 0.3"
        ),
        "equivariant_terms": [
            ("f", "r-1 = -0.7"),
            ("\\Delta f", "-2.0"),
            ("\\Delta^2 f", "-1.0"),
            ("f^3", "-1.0"),
        ],
        "standard_terms": [
            ("f", True, True, True, True),
            ("f^3", True, True, True, True),
            ("\\partial_{xx} f", True, True, True, True),
            ("\\partial_{yy} f", True, True, True, True),
            ("\\partial_{xxxx} f", True, True, True, True),
            ("\\partial_{xxyy} f", True, True, True, True),
            ("\\partial_{yyyy} f", True, True, True, True),
        ],
        "notes": (
            "Swift-Hohenberg contains a 4th-order biharmonic term "
            "$\\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ which "
            "is the main source of difficulty.  All four libraries "
            "can represent every term individually, but the equivariant "
            "basis encodes the isotropic combination as a single operator, "
            "giving it a $3\\times$ sparsity advantage over coordinate bases "
            "for the biharmonic term."
        ),
    },
    "cahn_hilliard": {
        "title": "Cahn-Hilliard equation",
        "latex": (
            "f_t = -\\Delta^2 f + \\Delta(f^3 - f) "
            "= -\\Delta f + 3f^2\\Delta f + 6f\\|\\nabla f\\|^2 "
            "- \\varepsilon^2\\Delta^2 f, \\qquad \\varepsilon = 0.2"
        ),
        "equivariant_terms": [
            ("\\Delta f", "-1.0"),
            ("f^2\\cdot\\Delta f", "+3.0"),
            ("f\\cdot\\|\\nabla f\\|^2", "+6.0"),
            ("\\Delta^2 f", "-\\varepsilon^2 = -0.04"),
        ],
        "standard_terms": [
            ("\\partial_{xx} f", True, True, True, True),
            ("\\partial_{yy} f", True, True, True, True),
            ("f^2\\cdot\\partial_{xx} f", True, False, True, True),
            ("f^2\\cdot\\partial_{yy} f", True, False, True, True),
            ("f\\cdot f_x^2", False, False, True, True),
            ("f\\cdot f_y^2", False, False, True, True),
            ("\\partial_{xxxx} f", True, True, True, True),
            ("\\partial_{xxyy} f", True, True, True, True),
            ("\\partial_{yyyy} f", True, True, True, True),
        ],
        "notes": (
            "Cahn-Hilliard combines second and fourth-order terms with "
            "nonlinear products.  PDE-FIND cannot represent $f\\cdot f_x^2$ "
            "(products of derivatives).  WSINDy cannot represent either "
            "$f^2\\cdot f_{xx}$ or $f\\cdot f_x^2$ individually, but the "
            "combination $3f^2\\Delta f + 6f\\|\\nabla f\\|^2$ equals "
            "$\\Delta(f^3)$, which IS representable as "
            "$(f^3)_{xx} + (f^3)_{yy}$ in the WSINDy basis."
        ),
    },
    "ks": {
        "title": "Kuramoto-Sivashinsky equation",
        "latex": (
            "f_t = -\\Delta f - \\Delta^2 f - "
            "\\tfrac12\\|\\nabla f\\|^2"
        ),
        "equivariant_terms": [
            ("\\Delta f", "-1.0"),
            ("\\Delta^2 f", "-1.0"),
            ("\\|\\nabla f\\|^2", "-0.5"),
        ],
        "standard_terms": [
            ("\\partial_{xx} f", True, True, True, True),
            ("\\partial_{yy} f", True, True, True, True),
            ("\\partial_{xxxx} f", True, True, True, True),
            ("\\partial_{xxyy} f", True, True, True, True),
            ("\\partial_{yyyy} f", True, True, True, True),
            ("f_x^2", False, False, True, True),
            ("f_y^2", False, False, True, True),
        ],
        "notes": (
            "KS contains the gradient norm $\\|\\nabla f\\|^2 = f_x^2 + f_y^2$ "
            "which neither PDE-FIND nor WSINDy can represent as a sparse "
            "column (PDE-FIND lacks products of derivatives; in WSINDy, "
            "$f_x^2$ is locked with $f\\cdot f_{xx}$ inside "
            "$\\partial_{xx}(f^2)$)."
        ),
    },
    "porous_medium": {
        "title": "Porous Medium equation",
        "latex": (
            "f_t = \\nabla\\cdot(f^2\\nabla f) = "
            "f^2\\Delta f + 2f\\|\\nabla f\\|^2"
        ),
        "equivariant_terms": [
            ("f^2\\cdot\\Delta f", "+1.0"),
            ("f\\cdot\\|\\nabla f\\|^2", "+2.0"),
        ],
        "standard_terms": [
            ("f^2\\cdot\\partial_{xx} f", True, True, True, True),
            ("f^2\\cdot\\partial_{yy} f", True, True, True, True),
            ("f\\cdot f_x^2", False, True, True, True),
            ("f\\cdot f_y^2", False, True, True, True),
        ],
        "notes": (
            "PDE-FIND cannot represent $f\\cdot f_x^2$ "
            "or $f\\cdot f_y^2$ (products of derivatives), so its "
            "minimum fn is 2."
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def best_in_group(group_rows):
    def key(r):
        exact = 0 if r["exact_recovery"] == "1" else 1
        fp = float(r["fp"]) if r["fp"] != "nan" else float("inf")
        fn = float(r["fn"]) if r["fn"] != "nan" else float("inf")
        e = float(r["E"]) if r["E"] != "nan" else float("inf")
        return (exact, fp, fn, e)
    return min(group_rows, key=key)


def fmt_val(v, d=4):
    try:
        return f"{float(v):.{d}f}"
    except (ValueError, TypeError):
        return str(v)


def is_failure(r):
    if r is None:
        return True
    return int(r["fn"]) >= 3 or float(r["E"]) >= 0.99


def status_mark(r):
    if r is None:
        return ""
    if r["exact_recovery"] == "1":
        return " \u2713"
    if is_failure(r):
        return " \u2717"
    return ""


def eq_wins(eq, ne):
    """Equivariant vs non-equivariant: 🟢 clear eq win, 🟩 marginal, ⚪ tie,
    🟥 marginal non-eq win, 🔴 clear non-eq win, — both failed."""
    if eq is None and ne is None:
        return ""
    eq_fail = is_failure(eq)
    ne_fail = is_failure(ne)
    if eq_fail and ne_fail:
        return "\u2014"
    if eq_fail:
        return "\U0001f534"
    if ne_fail:
        return "\U0001f7e2"
    eq_tot = int(eq["fp"]) + int(eq["fn"])
    ne_tot = int(ne["fp"]) + int(ne["fn"])
    if eq_tot < ne_tot:
        return "\U0001f7e2"
    if ne_tot < eq_tot:
        return "\U0001f534"
    eq_E = float(eq["E"])
    ne_E = float(ne["E"])
    if abs(eq_E - ne_E) < 1e-4:
        return "\u26aa"
    if eq_E < ne_E:
        return "\U0001f7e9"
    return "\U0001f7e5"


def get_meta(eq_name):
    for prefix, m in EQUATION_META.items():
        if eq_name.startswith(prefix):
            return m
    return None


# ═══════════════════════════════════════════════════════════════════════════
# Table builders
# ═══════════════════════════════════════════════════════════════════════════

LEGEND = """### Legend for summary table

| Symbol | Meaning |
|---|---|
| \U0001f7e2 | Equivariant **clear** win (fewer total fp+fn) |
| \U0001f7e9 | Equivariant **marginal** win (same fp+fn, lower E) |
| \u26aa | Tie (same fp+fn, same E) |
| \U0001f7e5 | Non-equivariant **marginal** win (same fp+fn, lower E) |
| \U0001f534 | Non-equivariant **clear** win (fewer total fp+fn) |
| \u2014 | Both methods failed, no comparison |
| \u2713 | Exact recovery (fp = fn = 0) |
| \u2717 | Failed |
"""


def best_eq(rows, noise, method, backend, profile):
    """Best equivariant result for a given config."""
    pt, pd = profile
    grp = [
        r for r in rows
        if r.get("profile_type", "") == pt
        and r.get("profile_direction", "") == pd
        and r["noise"] == noise and r["method"] == method
        and r["backend"] == backend and r["library"] == "equivariant"
        and r["fp"] != "nan"
    ]
    return best_in_group(grp) if grp else None


def best_ne(rows, noise, method, backend, profile):
    """Best non-equivariant result (excl. equivariant and mwsindy)."""
    pt, pd = profile
    grp = [
        r for r in rows
        if r.get("profile_type", "") == pt
        and r.get("profile_direction", "") == pd
        and r["noise"] == noise and r["method"] == method
        and r["backend"] == backend
        and r["library"] not in ("equivariant", "mwsindy")
        and r["fp"] != "nan"
    ]
    return best_in_group(grp) if grp else None


def best_mw(rows, noise, method, backend, profile):
    """Best MWSINDy result."""
    pt, pd = profile
    grp = [
        r for r in rows
        if r.get("profile_type", "") == pt
        and r.get("profile_direction", "") == pd
        and r["noise"] == noise and r["method"] == method
        and r["backend"] == backend and r["library"] == "mwsindy"
        and r["fp"] != "nan"
    ]
    return best_in_group(grp) if grp else None


def build_summary_table(rows):
    """Equivariant vs best non-equivariant, per anisotropic profile.

    Same format as build_summary.py but with a Profile column and only
    anisotropic noise profiles."""
    lines = [
        "## Summary: equivariant vs best non-equivariant (anisotropic noise)",
        "",
        "| Noise | Method | Backend | Profile | Eq best | Multigraphs fp/fn | Multigraphs E | Non-eq best | Non-eq fp/fn | Non-eq E | MWSINDy fp/fn | MWSINDy E |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]

    for noise, noise_label in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
        for method, method_label in METHOD_MAP.items():
            for be, be_label in [("spectral", "spectral"), ("polynomial", "polynomial")]:
                for pt, pd in PROFILE_ORDER:
                    label = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")
                    eq = best_eq(rows, noise, method, be, (pt, pd))
                    ne = best_ne(rows, noise, method, be, (pt, pd))
                    mw = best_mw(rows, noise, method, be, (pt, pd))
                    if eq is None and ne is None and mw is None:
                        continue

                    winner = eq_wins(eq, ne)

                    eq_fpf = f"{eq['fp']}/{eq['fn']}{status_mark(eq)}" if eq else "\u2014"
                    eq_E = fmt_val(eq["E"]) if eq else "\u2014"

                    if ne is None:
                        ne_lib = "\u2014"
                        ne_fpf = "\u2014"
                        ne_E = "\u2014"
                    else:
                        ne_lib = "\u2014" if is_failure(ne) else LIB_MAP.get(ne["library"], "")
                        ne_fpf = f"{ne['fp']}/{ne['fn']}{status_mark(ne)}"
                        ne_E = fmt_val(ne["E"])

                    mw_fpf = f"{mw['fp']}/{mw['fn']}{status_mark(mw)}" if mw else "\u2014"
                    mw_E = fmt_val(mw["E"]) if mw else "\u2014"

                    row = f"| {noise_label} | {method_label} | {be_label} | {label} | {winner} | {eq_fpf} | {eq_E} | {ne_lib} | {ne_fpf} | {ne_E} | {mw_fpf} | {mw_E} |"
                    lines.append(row)
    return "\n".join(lines)


def build_detailed_table(rows, method, method_label, backend, backend_label):
    """Detailed per-library table for one (method, backend)."""
    lines = [
        f"### {backend_label} backend",
        "",
        "| Library | Noise | Algorithm | d_tol/\u03b1 | Profile | fp/fn | n_active | E | Time (s) | Exact? |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]

    for lib, lib_label in [
        ("standard", "PDE-FIND"),
        ("wsindy", "WSINDy"),
        ("polynomial_jet", "PolyJet"),
        ("equivariant", "Multigraphs"),
        ("mwsindy", "MWSINDy"),
    ]:
        for n, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
            for alg, alg_label in [("strridge", "STRidge"), ("sr3", "SR3"), ("mstls", "MSTLS")]:
                for pt, pd in PROFILE_ORDER:
                    plabel = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")
                    g = [
                        r for r in rows
                        if r.get("profile_type", "") == pt
                        and r.get("profile_direction", "") == pd
                        and r["library"] == lib and r["noise"] == n
                        and r["method"] == method
                        and r["backend"] == backend and r["algorithm"] == alg
                        and r["fp"] != "nan"
                    ]
                    if not g:
                        continue
                    b = best_in_group(g)
                    fpf = f"{b['fp']}/{b['fn']}"
                    ex = "\U0001f7e2" if b["exact_recovery"] == "1" else "\U0001f534"
                    E = fmt_val(b["E"])
                    t = fmt_val(float(b.get("time_ms", 0)) / 1000, 1)
                    p = b.get("mstls_alpha", "") if alg == "mstls" else b.get("d_tol", "")
                    row = f"| {lib_label} | {nl} | {alg_label} | {p} | {plabel} | {fpf} | {b['n_active']} | {E} | {t} | {ex} |"
                    lines.append(row)
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════
# Main
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

        rep_lines = [
            "| Operator | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |",
            "|---|---|---|---|---|---|",
        ]
        for term, pf, ws, pj, mg in meta["standard_terms"]:
            def yn(b):
                return "\U0001f7e2" if b else "\U0001f534"
            rep_lines.append(
                f"| ${term}$ | | {yn(pf)} | {yn(ws)} | {yn(pj)} | {yn(mg)} |"
            )
        rep_table = "\n".join(rep_lines)

        eq_terms_lines = [
            "| Operator | True coefficient |",
            "|---|---|",
        ]
        for term, coeff in meta["equivariant_terms"]:
            eq_terms_lines.append(f"| ${term}$ | ${coeff}$ |")
        eq_terms_table = "\n".join(eq_terms_lines)

        methods_str = ", ".join(METHOD_MAP.get(m, m) for m in methods)
        algs_str = ", ".join(algorithms)
        libs_str = ", ".join(LIB_MAP.get(l, l) for l in libraries)
        profiles_str = ", ".join(
            PROFILE_LABELS.get(p, f"{p[0]}/{p[1]}") for p in profiles_found
        )

        n_libs = len(libraries)
        lib_plural = "ies" if n_libs > 1 else "y"

        md = f"""# {meta["title"]} \u2014 Anisotropic noise

The equation under consideration is

$${meta["latex"]}$$

which can be written in terms of $O(2)$-equivariant terms

{eq_terms_table}

or standard terms (PolyJet basis)

{rep_table}

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

{build_summary_table(eq_rows)}

## Pointwise formulation

{build_detailed_table(eq_rows, 'pointwise', 'pointwise', 'spectral', 'Spectral derivatives')}

{build_detailed_table(eq_rows, 'pointwise', 'pointwise', 'polynomial', 'Polynomial stencils')}

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


def _read_log_footer(csv_path):
    log_path = csv_path.replace(".csv", ".log")
    if not os.path.exists(log_path):
        return ""
    info = {}
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or ":" not in line:
                continue
            key, val = line.split(":", 1)
            info[key.strip()] = val.strip()
    if not info:
        return ""
    total = info.get("total_wall_time_min", "?")
    cores = info.get("omp_num_threads", "?")
    cpu = info.get("cpu_model", "?")
    return f"\n*Run: {total} min ({cpu}), OMP_NUM_THREADS={cores}*"


def main():
    parser = argparse.ArgumentParser(
        description="Build anisotropic-noise summary markdown from a sweep CSV."
    )
    parser.add_argument("csv", help="Path to sweep CSV file.")
    parser.add_argument(
        "--equation", "-e", default=None,
        help="Restrict to one equation (default: all in CSV).",
    )
    args = parser.parse_args()
    build_markdown(args.csv, args.equation)


if __name__ == "__main__":
    main()
