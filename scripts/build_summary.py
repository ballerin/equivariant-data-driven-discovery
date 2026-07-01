#!/usr/bin/env python3
"""
Build a per-equation summary markdown file from a sweep CSV.

Reads the unified sweep output (produced by scripts/sweep.py), selects the
best hyperparameters for each (library, noise, method, backend, algorithm)
group, and writes a markdown file with:

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
    python scripts/build_summary.py results/sweep.csv --equation porous_medium

The output is written to summaries/<equation>.md.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

# ── paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SUMMARIES_DIR = os.path.join(CODE_DIR, "summaries")
os.makedirs(SUMMARIES_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Per-equation configuration
# ═══════════════════════════════════════════════════════════════════════════

# Map library keys to display names
LIB_MAP = {
    "equivariant": "Multigraphs",
    "standard": "PDE-FIND",
    "polynomial_jet": "PolyJet",
    "wsindy": "WSINDy",
    "wsindy_proj": "WSINDy → O(2)",
    "mwsindy": "MWSINDy",
}

# Equation metadata: (title, latex, equivariant_terms, standard_terms, notes)
EQUATION_META = {
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
            # (term, pde_find, wsindy, poly_jet, multigraph)
            ("f^2\\cdot\\partial_{xx} f", True, True, True, True),
            ("f^2\\cdot\\partial_{yy} f", True, True, True, True),
            ("f\\cdot f_x^2", False, True, True, True),
            ("f\\cdot f_y^2", False, True, True, True),
        ],
        "notes": (
            "PDE-FIND cannot represent $f\\cdot f_x^2$ "
            "or $f\\cdot f_y^2$ (products of derivatives), so its "
            "minimum fn is 2.  "
            "WSINDy represents PME exactly via the chain rule: "
            "$f_t = \\frac13\\partial_x^2(f^3) + "
            "\\frac13\\partial_y^2(f^3)$, which "
            "expands to all 4 jet terms."
        ),
    },
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
            # (term, pde_find, wsindy, poly_jet, multigraph)
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
            # (term, pde_find, wsindy, poly_jet, multigraph)
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
            "giving it a 3\\times sparsity advantage over coordinate bases "
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
            "nonlinear products.  PDE-FIND can represent $f^2\\cdot f_{xx}$ "
            "and the pure derivatives but cannot represent $f\\cdot f_x^2$ "
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
            "$\\partial_{xx}(f^2)$).  Three KS regimes are simulated: "
            "cellular ($L=4\\pi$), weak-chaotic ($L=8\\pi$), "
            "chaotic ($L=16\\pi$)."
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def best_in_group(group_rows):
    """Lexicographic best: exact → fp → fn → E."""
    def key(r):
        exact = 0 if r["exact_recovery"] == "1" else 1
        fp = float(r["fp"]) if r["fp"] != "nan" else float("inf")
        fn = float(r["fn"]) if r["fn"] != "nan" else float("inf")
        e = float(r["E"]) if r["E"] != "nan" else float("inf")
        return (exact, fp, fn, e)
    return min(group_rows, key=key)


def fmt(v, d=4):
    """Format a float or return the string as-is."""
    try:
        return f"{float(v):.{d}f}"
    except (ValueError, TypeError):
        return str(v)


def is_failure(r):
    """A result is a failure if it loses enough true terms (fn >= 3)
    or the coefficient error is near 100% (E >= 0.99)."""
    if r is None:
        return True
    fn = int(r["fn"])
    E = float(r["E"])
    return fn >= 3 or E >= 0.99


def is_exact(r):
    """Check for exact recovery (fp=0, fn=0)."""
    return r is not None and r["exact_recovery"] == "1"


def status_mark(r):
    """Return subtle symbol: ✓ exact, ✗ failure, '' otherwise."""
    if r is None:
        return ""
    if r["exact_recovery"] == "1":
        return " \u2713"
    if is_failure(r):
        return " \u2717"
    return ""


def best_eq(rows, noise, method, backend):
    """Best equivariant (direct multigraph) result."""
    grp = [
        r for r in rows
        if r["noise"] == noise and r["method"] == method
        and r["backend"] == backend and r["library"] == "equivariant"
        and r["fp"] != "nan"
    ]
    return best_in_group(grp) if grp else None


def best_mwsindy(rows, noise, method, backend):
    """Best MWSINDy result (WSINDy regression + multigraph fine-tuning)."""
    grp = [
        r for r in rows
        if r["noise"] == noise and r["method"] == method
        and r["backend"] == backend and r["library"] == "mwsindy"
        and r["fp"] != "nan"
    ]
    return best_in_group(grp) if grp else None


def best_non_eq(rows, noise, method, backend):
    """Best non-equivariant result (excludes Multigraph and MWSINDy)."""
    grp = [
        r for r in rows
        if r["noise"] == noise and r["method"] == method
        and r["backend"] == backend
        and r["library"] not in ("equivariant", "mwsindy")
        and r["fp"] != "nan"
    ]
    return best_in_group(grp) if grp else None


def _eq_best(eq, mw):
    """Which equivariant method is better: direct multigraph or MWSINDy?

    Clear win = one method has better fp/fn (fewer total false positives +
    false negatives).  Marginal win = same fp/fn, lower E.  Tie = same fp/fn
    and same E.

    Returns: 🟢 MWSINDy clear win, 🟩 MWSINDy marginal,
             ⚪ tie / both fail, 🟥 Multigraph marginal, 🔴 Multigraph clear.
    """
    if eq is None and mw is None:
        return "\u2014"

    eq_fail = is_failure(eq)
    mw_fail = is_failure(mw)
    if eq_fail and mw_fail:
        return "\u26aa"
    if eq_fail:
        return "\U0001f7e2"   # MWSINDy clear win (eq failed)
    if mw_fail:
        return "\U0001f534"   # Multigraph clear win (mw failed)

    eq_total = int(eq["fp"]) + int(eq["fn"])
    mw_total = int(mw["fp"]) + int(mw["fn"])
    eq_E = float(eq["E"])
    mw_E = float(mw["E"])

    if eq_total < mw_total:
        return "\U0001f534"   # Multigraph clear win
    elif mw_total < eq_total:
        return "\U0001f7e2"   # MWSINDy clear win
    else:
        # Same fp/fn — decide by E (marginal if close, tie if identical)
        if abs(eq_E - mw_E) < 1e-4:
            return "\u26aa"
        elif mw_E < eq_E:
            return "\U0001f7e9"   # MWSINDy marginal
        else:
            return "\U0001f7e5"   # Multigraph marginal


def eq_wins(eq, ne):
    """5-level winner indicator: equivariant vs best non-equivariant.

    Clear win  = one method has strictly fewer total fp+fn than the other.
    Marginal   = same fp+fn totals, but one method has lower E.
    Tie        = same fp+fn and same E (within 1e-4).

    Returns
    -------
    str
        🟢 clear eq win, 🟩 marginal eq win, ⚪ tie,
        🟥 marginal non-eq win, 🔴 clear non-eq win, — both failed.
    """
    if eq is None and ne is None:
        return ""

    eq_fail = is_failure(eq)
    ne_fail = is_failure(ne)

    if eq_fail and ne_fail:
        return "\u2014"
    if eq_fail:
        return "\U0001f534"   # non-eq clear win
    if ne_fail:
        return "\U0001f7e2"   # eq clear win

    eq_total = int(eq["fp"]) + int(eq["fn"])
    ne_total = int(ne["fp"]) + int(ne["fn"])
    eq_E = float(eq["E"])
    ne_E = float(ne["E"])

    if eq_total < ne_total:
        return "\U0001f7e2"   # eq clear win
    elif ne_total < eq_total:
        return "\U0001f534"   # non-eq clear win
    else:
        # Same fp/fn — decide by E (marginal)
        if abs(eq_E - ne_E) < 1e-4:
            return "\u26aa"
        elif eq_E < ne_E:
            return "\U0001f7e9"   # eq marginal
        else:
            return "\U0001f7e5"   # non-eq marginal


# ═══════════════════════════════════════════════════════════════════════════
# Table builders
# ═══════════════════════════════════════════════════════════════════════════

def build_summary_table(rows, eq_meta):
    """Build the equivariant-vs-best-non-equivariant summary table."""
    lines = [
        "## Summary: equivariant vs best non-equivariant",
        "",
        "| Noise | Method | Backend | Eq best | Multigraphs fp/fn | "
        "Multigraphs E | Non-eq best | Non-eq fp/fn | Non-eq E | "
        "MWSINDy fp/fn | MWSINDy E |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for noise, noise_label in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
        for method, method_label in [
            ("pointwise", "pointwise"),
            ("weak_direct", "weak (no IBP)"),
            ("weak_ibp", "weak (IBP)"),
        ]:
            for be, be_label in [("spectral", "spectral"), ("polynomial", "polynomial")]:
                eq = best_eq(rows, noise, method, be)
                mw = best_mwsindy(rows, noise, method, be)
                ne = best_non_eq(rows, noise, method, be)
                if eq is None and mw is None and ne is None:
                    continue

                # Winner: Multigraphs vs best non-equivariant (MWSINDy is
                # reported but not used in the comparison).
                winner = eq_wins(eq, ne)

                eq_fpf = f"{eq['fp']}/{eq['fn']}{status_mark(eq)}" if eq else "\u2014"
                eq_E = fmt(eq["E"]) if eq else "\u2014"

                # Non-eq: if no result at all, blank.  If a result exists
                # but is a failure, show the library as — but still
                # report the fp/fn and E so the reader can see what
                # happened.
                if ne is None:
                    ne_lib = "\u2014"
                    ne_fpf = "\u2014"
                    ne_E = "\u2014"
                else:
                    ne_lib = "\u2014" if is_failure(ne) else LIB_MAP.get(ne["library"], "")
                    ne_fpf = f"{ne['fp']}/{ne['fn']}{status_mark(ne)}"
                    ne_E = fmt(ne["E"])

                # MWSINDy — data columns only, not part of the comparison
                mw_fpf = f"{mw['fp']}/{mw['fn']}{status_mark(mw)}" if mw else "\u2014"
                mw_E = fmt(mw["E"]) if mw else "\u2014"

                lines.append(
                    f"| {noise_label} | {method_label} | {be_label} | "
                    f"{winner} | {eq_fpf} | {eq_E} | "
                    f"{ne_lib} | {ne_fpf} | {ne_E} | "
                    f"{mw_fpf} | {mw_E} |"
                )
    return "\n".join(lines)


def build_pointwise_table(rows, backend, backend_label):
    """Detailed pointwise table for one backend."""
    lines = [
        f"### {backend_label} backend",
        "",
        "| Library | Noise | Algorithm | d_tol/\u03b1 | "
        "fp/fn | n_active | E | Time (s) | Exact? |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for lib, lib_label in [
        ("standard", "PDE-FIND"),
        ("wsindy", "WSINDy"),
        ("polynomial_jet", "PolyJet"),
        ("equivariant", "Multigraphs"),
        ("mwsindy", "MWSINDy"),
    ]:
        for n, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
            for a, al in [("strridge", "STRidge"), ("sr3", "SR3"), ("mstls", "MSTLS")]:
                g = [
                    r for r in rows
                    if r["library"] == lib and r["noise"] == n
                    and r["method"] == "pointwise"
                    and r["backend"] == backend and r["algorithm"] == a
                    and r["fp"] != "nan"
                ]
                if not g:
                    continue
                b = best_in_group(g)
                fpf = f"{b['fp']}/{b['fn']}"
                ex = "\U0001f7e2" if b["exact_recovery"] == "1" else "\U0001f534"
                E = fmt(b["E"])
                t = fmt(float(b.get("time_ms", 0)) / 1000, 1)
                p = b.get("mstls_alpha", "") if a == "mstls" else b.get("d_tol", "")
                lines.append(
                    f"| {lib_label} | {nl} | {al} | {p} | "
                    f"{fpf} | {b['n_active']} | {E} | {t} | {ex} |"
                )
    return "\n".join(lines)


def build_weak_table(rows, backend, backend_label):
    """Detailed weak-form table for one backend."""
    lines = [
        f"### {backend_label} backend",
        "",
        "| Library | Noise | Algorithm | d_tol/\u03b1 | m | IBP | "
        "fp/fn | n_active | E | Time (s) | Exact? |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for lib, lib_label in [
        ("standard", "PDE-FIND"),
        ("wsindy", "WSINDy"),
        ("polynomial_jet", "PolyJet"),
        ("equivariant", "Multigraphs"),
        ("mwsindy", "MWSINDy"),
    ]:
        for n, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
            for mtd, ibp in [("weak_direct", "no"), ("weak_ibp", "yes")]:
                for a, al in [("strridge", "STRidge"), ("sr3", "SR3"), ("mstls", "MSTLS")]:
                    g = [
                        r for r in rows
                        if r["library"] == lib and r["noise"] == n
                        and r["method"] == mtd
                        and r["backend"] == backend and r["algorithm"] == a
                        and r["fp"] != "nan"
                    ]
                    if not g:
                        continue
                    b = best_in_group(g)
                    fpf = f"{b['fp']}/{b['fn']}"
                    ex = "\U0001f7e2" if b["exact_recovery"] == "1" else "\U0001f534"
                    E = fmt(b["E"])
                    t = fmt(float(b.get("time_ms", 0)) / 1000, 1)
                    p = b.get("mstls_alpha", "") if a == "mstls" else b.get("d_tol", "")
                    lines.append(
                        f"| {lib_label} | {nl} | {al} | {p} | "
                        f"{b.get('weak_m', '')} | {ibp} | "
                        f"{fpf} | {b['n_active']} | {E} | {t} | {ex} |"
                    )
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════
# Legend
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


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def build_markdown(csv_path, equation=None):
    """Read a sweep CSV and return the complete markdown string."""
    with open(csv_path) as f:
        all_rows = list(csv.DictReader(f))

    # Filter to one equation if requested
    eq_names = set(r["equation"] for r in all_rows)
    if equation:
        if equation not in eq_names:
            # Try prefix match
            matches = [e for e in eq_names if e.startswith(equation)]
            if not matches:
                print(f"Equation {equation!r} not found. Available: {sorted(eq_names)}")
                sys.exit(1)
            equation = matches[0]
        rows = [r for r in all_rows if r["equation"] == equation]
    else:
        # Use the first equation found
        equation = sorted(eq_names)[0]
        rows = [r for r in all_rows if r["equation"] == equation]

    # Look up metadata
    meta = None
    for prefix, m in EQUATION_META.items():
        if equation.startswith(prefix):
            meta = m
            break
    if meta is None:
        print(f"No metadata for equation {equation!r}. "
              f"Available: {list(EQUATION_META)}")
        sys.exit(1)

    # Build representability table
    rep_lines = [
        f"| Operator | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |",
        "|---|---|---|---|---|---|",
    ]
    for term, pf, ws, pj, mg in meta["standard_terms"]:
        def yn(b):
            return "\U0001f7e2" if b else "\U0001f534"
        rep_lines.append(
            f"| ${term}$ | | {yn(pf)} | {yn(ws)} | {yn(pj)} | {yn(mg)} |"
        )
    rep_table = "\n".join(rep_lines)

    # Build equivariant terms table
    eq_terms_lines = [
        "| Operator | True coefficient |",
        "|---|---|",
    ]
    for term, coeff in meta["equivariant_terms"]:
        eq_terms_lines.append(f"| ${term}$ | ${coeff}$ |")
    eq_terms_table = "\n".join(eq_terms_lines)

    # Get sweep file info
    sweep_name = os.path.basename(csv_path)

    # Assemble full markdown
    md = f"""# {meta["title"]}

The equation under consideration is

$${meta["latex"]}$$

which can be written in terms of $O(2)$-equivariant terms

{eq_terms_table}

or standard terms (PolyJet basis)

{rep_table}

All libraries are evaluated in the common polynomial-jet basis (105 monomials,
max order 4, max degree 3).  fp/fn/E are computed against the same 4-term
ground truth for every library.  {meta["notes"]}

Sweep configuration: all 4 libraries, 3 algorithms (STRidge, SR3, MSTLS),
3 methods (pointwise, weak-direct, weak-IBP), 2 backends (spectral, polynomial),
3 noise levels (clean, 2%, 5%).

{LEGEND}

{build_summary_table(rows, meta)}

## Pointwise formulation

{build_pointwise_table(rows, 'spectral', 'Spectral derivatives')}

{build_pointwise_table(rows, 'polynomial', 'Polynomial stencils')}

## Weak formulation

{build_weak_table(rows, 'spectral', 'Spectral derivatives')}

{build_weak_table(rows, 'polynomial', 'Polynomial stencils')}

---
*Source: `{sweep_name}`*
{_read_log_footer(csv_path)}
"""
    return md


def _read_log_footer(csv_path):
    """Read companion .log file and return a footer string with timing info."""
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
    return (
        f"\n*Run: {total} min ({cpu}), "
        f"OMP_NUM_THREADS={cores}*"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Build per-equation summary markdown from a sweep CSV."
    )
    parser.add_argument(
        "csv", help="Path to sweep CSV file (e.g. results/sweep_20260619_170245_pme_v3.csv)"
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
        with open(args.csv) as f:
            reader = csv.DictReader(f)
            eq_names = set(r["equation"] for r in reader)
        eq = args.equation or sorted(eq_names)[0]
        out_path = os.path.join(SUMMARIES_DIR, f"{eq}.md")

    with open(out_path, "w") as f:
        f.write(md)

    print(f"Written {out_path}")


if __name__ == "__main__":
    main()
