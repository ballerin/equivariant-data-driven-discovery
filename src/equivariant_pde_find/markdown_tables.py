"""
Programmatic markdown table builder and shared helper functions for
sweep summarization scripts.

Provides
--------
MarkdownTable
    A safe builder for pipe-delimited markdown tables.  Rows are added
    one at a time; the render method handles alignment and escaping.
best_in_group
    Lexicographic best-row selection from a group of sweep rows.
Shared helpers
    fmt_val, is_failure, status_mark, eq_wins, _read_log_footer.
"""

import math
import os
from collections import defaultdict


# ═══════════════════════════════════════════════════════════════════════════
# MarkdownTable — safe, programmatic table builder
# ═══════════════════════════════════════════════════════════════════════════

class MarkdownTable:
    """Builder for pipe-delimited markdown tables.

    Use this instead of manual f-string concatenation to avoid malformed
    output from unescaped pipes, missing delimiters, or misaligned columns.

    Parameters
    ----------
    headers : list of str
        Column header names.
    align : list of str, optional
        One of ``"left"``, ``"center"``, ``"right"`` per column.
        Default is left-aligned for all columns.

    Example
    -------
    >>> t = MarkdownTable(["Library", "fp/fn", "E"], align=["left", "center", "right"])
    >>> t.add_row("Multigraphs", "0/0", "0.0012")
    >>> t.add_row("PDE-FIND", "3/1", "0.3400")
    >>> print(t.render())
    | Library      | fp/fn |      E |
    |--------------|:-----:|-------:|
    | Multigraphs  |  0/0  | 0.0012 |
    | PDE-FIND     |  3/1  | 0.3400 |
    """

    def __init__(self, headers, align=None):
        self.headers = list(headers)
        self._rows = []
        if align is None:
            self._align = ["left"] * len(headers)
        else:
            if len(align) != len(headers):
                raise ValueError(
                    f"align has {len(align)} entries but headers has "
                    f"{len(headers)}"
                )
            for a in align:
                if a not in ("left", "center", "right"):
                    raise ValueError(
                        f"Invalid alignment {a!r}; use 'left', 'center', "
                        f"or 'right'"
                    )
            self._align = list(align)

    def add_row(self, *cells):
        """Add a data row.  Each cell is converted to str."""
        if len(cells) != len(self.headers):
            raise ValueError(
                f"Expected {len(self.headers)} cells, got {len(cells)}"
            )
        self._rows.append([str(c) for c in cells])

    # ── rendering ──────────────────────────────────────────────────────

    def render(self):
        """Render the full table as a multi-line string."""
        all_rows = [self.headers] + self._rows
        col_widths = self._compute_widths(all_rows)

        lines = []
        lines.append(self._render_row(self.headers, col_widths))
        lines.append(self._render_separator(col_widths))
        for row in self._rows:
            lines.append(self._render_row(row, col_widths))
        return "\n".join(lines)

    # ── internal helpers ───────────────────────────────────────────────

    @staticmethod
    def _compute_widths(rows):
        n_cols = len(rows[0])
        widths = [0] * n_cols
        for row in rows:
            for j, cell in enumerate(row):
                widths[j] = max(widths[j], len(cell))
        return widths

    def _render_separator(self, widths):
        parts = []
        for j, w in enumerate(widths):
            a = self._align[j]
            if a == "left":
                parts.append(":" + "-" * (w + 1))
            elif a == "center":
                parts.append(":" + "-" * w + ":")
            else:  # right
                parts.append("-" * (w + 1) + ":")
        return "|" + "|".join(parts) + "|"

    @staticmethod
    def _render_row(cells, widths):
        parts = []
        for j, (cell, w) in enumerate(zip(cells, widths)):
            parts.append(" " + cell.ljust(w) + " ")
        return "|" + "|".join(parts) + "|"


# ═══════════════════════════════════════════════════════════════════════════
# Shared helper functions
# ═══════════════════════════════════════════════════════════════════════════

def best_in_group(group_rows):
    """Lexicographic best: exact → fp → fn → E_coeff.

    Parameters
    ----------
    group_rows : list of dict
        Sweep rows with keys ``exact_recovery``, ``fp``, ``fn``, ``E``.

    Returns
    -------
    dict
        The best row.
    """
    def key(r):
        exact = 0 if str(r.get("exact_recovery", "0")) == "1" else 1
        fp = float(r["fp"]) if r.get("fp", "nan") != "nan" else float("inf")
        fn = float(r["fn"]) if r.get("fn", "nan") != "nan" else float("inf")
        e = float(r["E_coeff"]) if r.get("E_coeff", "nan") != "nan" else float("inf")
        return (exact, fp, fn, e)
    return min(group_rows, key=key)


def fmt_val(v, d=4):
    """Format a float to *d* decimal places, or return the string as-is."""
    try:
        return f"{float(v):.{d}f}"
    except (ValueError, TypeError):
        return str(v)


def is_failure(r):
    """A result is a failure if fp > 0 or fn > 0 (not exact recovery)."""
    if r is None:
        return True
    try:
        fp = int(r["fp"])
        fn = int(r["fn"])
    except (ValueError, TypeError, KeyError):
        return True
    return fp > 0 or fn > 0


def status_mark(r):
    """Return subtle symbol: ✓ exact, ✗ failure, '' otherwise."""
    if r is None:
        return ""
    if str(r.get("exact_recovery", "0")) == "1":
        return " \u2713"
    else:
        return " \u2717"
    return ""


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
    eq_E = float(eq["E_coeff"])
    ne_E = float(ne["E_coeff"])

    if eq_total < ne_total:
        return "\U0001f7e2"   # eq clear win
    elif ne_total < eq_total:
        return "\U0001f534"   # non-eq clear win
    else:
        if abs(eq_E - ne_E) < 1e-4:
            return "\u26aa"
        elif eq_E < ne_E:
            return "\U0001f7e9"   # eq marginal
        else:
            return "\U0001f7e5"   # non-eq marginal


def compare_bases(o, m):
    """Compare original vs modified multigraph (for multigraph comparison).

    Returns color indicator: 🟢 modified wins, 🔴 original wins, etc.
    """
    eo = int(o.get("exact_recovery", "0"))
    em = int(m.get("exact_recovery", "0"))
    if em and not eo:
        return "\U0001f7e2"
    if eo and not em:
        return "\U0001f534"
    so = int(o["fp"]) + int(o["fn"])
    sm = int(m["fp"]) + int(m["fn"])
    if sm < so:
        return "\U0001f7e2"
    if so < sm:
        return "\U0001f534"
    if sm == so and float(m["E_coeff"]) < float(o["E_coeff"]):
        return "\U0001f7e9"
    if sm == so and float(o["E_coeff"]) < float(m["E_coeff"]):
        return "\U0001f7e5"
    return "\u26aa"


def _read_log_info(csv_path):
    """Read companion .log file and return a dict of timing info.

    Returns empty dict if no log file exists or it's unparseable.
    """
    log_path = csv_path.replace(".csv", ".log")
    if not os.path.exists(log_path):
        return {}
    info = {}
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or ":" not in line:
                continue
            key, val = line.split(":", 1)
            info[key.strip()] = val.strip()
    return info


def _read_log_footer(csv_path):
    """Read companion .log file and return a footer string with timing info."""
    info = _read_log_info(csv_path)
    if not info:
        return ""
    total = info.get("total_wall_time_min", "?")
    cores = info.get("omp_num_threads", "?")
    cpu = info.get("cpu_model", "?")
    return (
        f"\n*Run: {total} min ({cpu}), "
        f"OMP_NUM_THREADS={cores}*"
    )


# ═══════════════════════════════════════════════════════════════════════════
# Higher-level table builders (used by build_summary, build_anisotropic_summary)
# ═══════════════════════════════════════════════════════════════════════════

def _best_in_lib(rows, lib, noise, method, backend, algorithm=None):
    """Filter rows by library/noise/method/backend/algorithm, return best."""
    grp = [
        r for r in rows
        if r["library"] == lib
        and str(r.get("noise", "")) == str(noise)
        and r["method"] == method
        and r["backend"] == backend
        and r.get("fp", "nan") != "nan"
    ]
    if algorithm:
        grp = [r for r in grp if r.get("algorithm", "") == algorithm]
    return best_in_group(grp) if grp else None


def _best_eq(rows, noise, method, backend):
    """Best multigraph result (excludes oracle)."""
    grp = [
        r for r in rows
        if r["library"] == "multigraphs"
        and r.get("algorithm", "") != "oracle"
        and str(r.get("noise", "")) == str(noise)
        and r["method"] == method
        and r["backend"] == backend
        and r.get("fp", "nan") != "nan"
    ]
    return best_in_group(grp) if grp else None


def _best_mwsindy(rows, noise, method, backend):
    """Best MWSINDy result."""
    return _best_in_lib(rows, "mwsindy", noise, method, backend)


def _best_non_eq(rows, noise, method, backend):
    """Best non-multigraph result (excludes multigraphs, mwsindy, oracle)."""
    grp = [
        r for r in rows
        if r["library"] not in ("multigraphs", "mwsindy")
        and r.get("algorithm", "") != "oracle"
        and str(r.get("noise", "")) == str(noise)
        and r["method"] == method
        and r["backend"] == backend
        and r.get("fp", "nan") != "nan"
    ]
    return best_in_group(grp) if grp else None


def _scalable_rows(rows):
    """Return rows excluding non-scalable algorithms (brute-force)."""
    return [r for r in rows if r.get("algorithm", "") != "brute"]


def _oracle_E(rows, noise, method, backend):
    """Best oracle E_data for a given (noise, method, backend).

    Returns the minimum E_data among multigraph oracle rows for this
    configuration — the irreducible data-fit floor for that method.
    """
    oracle_rows = [
        r for r in rows
        if r.get("algorithm") == "oracle"
        and r["library"] == "multigraphs"
        and str(r.get("noise", "")) == str(noise)
        and r["method"] == method
        and r["backend"] == backend
        and r.get("E_data", "") not in ("", "nan")
    ]
    if not oracle_rows:
        return None
    return min(float(r["E_data"]) for r in oracle_rows)


def build_standard_summary_table(rows):
    """Build the equivariant-vs-best-non-equivariant summary table.

    Only scalable algorithms (STRidge, MSTLS, RFE) are considered;
    brute-force exhaustive search is excluded from the summary.
    Adds an Oracle E column showing the best achievable residual
    (true coefficients plugged into the multigraph library).
    """
    from .summary_metadata import LIB_NAMES, ALGORITHM_LABELS

    rows = _scalable_rows(rows)

    t = MarkdownTable(
        [
            "Noise", "Formulation", "Backend", "Eq wins?",
            "Multigraphs fp/fn", "Mg E_coeff", "Mg E_data",
            "Std best", "Std fp/fn", "Std E_coeff", "Std E_data",
            "E_data floor",
        ],
        align=["left"] + ["center"] * 3 + ["left"] * 7 + ["right"],
    )

    for noise, noise_label in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
        for method, method_label in [
            ("pointwise", "pointwise"),
            ("weak_direct", "weak (no IBP)"),
            ("weak_ibp", "weak (IBP)"),
        ]:
            for be, be_label in [("spectral", "spectral"), ("polynomial", "polynomial")]:
                eq = _best_eq(rows, noise, method, be)
                ne = _best_non_eq(rows, noise, method, be)
                if eq is None and ne is None:
                    continue

                winner = eq_wins(eq, ne)

                eq_fpf = f"{eq['fp']}/{eq['fn']}{status_mark(eq)}" if eq else "\u2014"
                eq_E = fmt_val(eq["E_coeff"]) if eq else "\u2014"
                eq_E_data = fmt_val(eq.get("E_data", "")) if eq else "\u2014"

                if ne is None:
                    ne_lib = "\u2014"
                    ne_fpf = "\u2014"
                    ne_E = "\u2014"
                    ne_E_data = "\u2014"
                else:
                    ne_lib = "\u2014" if is_failure(ne) else (
                        f"{LIB_NAMES.get(ne['library'], '')}/"
                        f"{ALGORITHM_LABELS.get(ne.get('algorithm', ''), ne.get('algorithm', ''))}"
                    )
                    ne_fpf = f"{ne['fp']}/{ne['fn']}{status_mark(ne)}"
                    ne_E = fmt_val(ne["E_coeff"])
                    ne_E_data = fmt_val(ne.get("E_data", ""))

                oracle_e = _oracle_E(rows, noise, method, be)
                oracle_str = fmt_val(oracle_e) if oracle_e is not None else "\u2014"

                t.add_row(
                    noise_label, method_label, be_label, winner,
                    eq_fpf, eq_E, eq_E_data, ne_lib, ne_fpf, ne_E, ne_E_data, oracle_str,
                )
    return t


def build_pointwise_table(rows, backend, backend_label):
    """Detailed pointwise table for one backend."""
    t = MarkdownTable(
        [
            "Library", "Noise", "Algorithm", "d_tol/\u03b1",
            "fp/fn", "n_active", "E_coeff", "E_data", "Time (s)", "Exact?",
        ],
    )

    for lib, lib_label in [
        ("standard", "PDE-FIND"),
        ("wsindy", "WSINDy"),
        ("polynomial_jet", "PolyJet"),
        ("multigraphs", "Multigraphs"),
    ]:
        for n, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
            for a, al in [("strridge", "STRidge"), ("mstls", "MSTLS"), ("rfe", "RFE"), ("brute", "Brute-force")]:
                b = _best_in_lib(rows, lib, n, "pointwise", backend, a)
                if b is None:
                    continue
                fpf = f"{b['fp']}/{b['fn']}"
                ex = "\U0001f7e2" if str(b.get("exact_recovery", "0")) == "1" else "\U0001f534"
                E = fmt_val(b["E_coeff"])
                E_data = fmt_val(b.get("E_data", ""))
                tim = fmt_val(float(b.get("time_ms", 0)) / 1000, 1)
                if a == "mstls":
                    p = b.get("mstls_alpha", "")
                elif a == "brute":
                    p = f"L₀={b.get('d_tol', '')}"
                else:
                    p = b.get("d_tol", "")
                t.add_row(lib_label, nl, al, str(p), fpf, str(b["n_active"]), E, E_data, tim, ex)
    return t


def build_weak_table(rows, backend, backend_label):
    """Detailed weak-form table for one backend."""
    t = MarkdownTable(
        [
            "Library", "Noise", "Algorithm", "d_tol/\u03b1", "m", "IBP",
            "fp/fn", "n_active", "E_coeff", "E_data", "Time (s)", "Exact?",
        ],
    )

    for lib, lib_label in [
        ("standard", "PDE-FIND"),
        ("wsindy", "WSINDy"),
        ("polynomial_jet", "PolyJet"),
        ("multigraphs", "Multigraphs"),
    ]:
        for n, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
            for mtd, ibp in [("weak_direct", "no"), ("weak_ibp", "yes")]:
                for a, al in [("strridge", "STRidge"), ("mstls", "MSTLS"), ("rfe", "RFE"), ("brute", "Brute-force")]:
                    b = _best_in_lib(rows, lib, n, mtd, backend, a)
                    if b is None:
                        continue
                    fpf = f"{b['fp']}/{b['fn']}"
                    ex = "\U0001f7e2" if str(b.get("exact_recovery", "0")) == "1" else "\U0001f534"
                    E = fmt_val(b["E_coeff"])
                    E_data = fmt_val(b.get("E_data", ""))
                    tim = fmt_val(float(b.get("time_ms", 0)) / 1000, 1)
                    if a == "mstls":
                        p = b.get("mstls_alpha", "")
                    elif a == "brute":
                        p = f"L₀={b.get('d_tol', '')}"
                    else:
                        p = b.get("d_tol", "")
                    t.add_row(
                        lib_label, nl, al, str(p), str(b.get("weak_m", "")),
                        ibp, fpf, str(b["n_active"]), E, E_data, tim, ex,
                    )
    return t


def build_oracle_table(rows):
    """Build a table of oracle E_data per configuration.

    The oracle plugs true multigraph coefficients into each method's
    (Theta, Ft) pair, giving per-configuration noise floors.
    """
    t = MarkdownTable(
        ["Noise", "Formulation", "Backend", "m", "IBP", "E_data floor", "cond(Θ)"],
        align=["left", "left", "left", "center", "center", "right", "left"],
    )

    oracle_rows = [
        r for r in rows
        if r.get("algorithm") == "oracle"
        and r["library"] == "multigraphs"
        and r.get("E_data", "") not in ("", "nan")
    ]
    if not oracle_rows:
        return None

    for noise, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
        for method, ml in [("pointwise", "pointwise"),
                           ("weak_direct", "weak (no IBP)"),
                           ("weak_ibp", "weak (IBP)")]:
            for be, bl in [("spectral", "spectral"), ("polynomial", "polynomial")]:
                if method == "pointwise":
                    m_vals = [0]
                else:
                    m_vals = [4, 8, 12]
                for mv in m_vals:
                    ibp = "yes" if method == "weak_ibp" else (
                        "no" if method == "weak_direct" else "")
                    matches = [
                        r for r in oracle_rows
                        if str(r["noise"]) == str(noise)
                        and r["method"] == method
                        and r["backend"] == be
                        and str(r.get("weak_m", "0")) == str(mv)
                    ]
                    if not matches:
                        continue
                    best = min(matches, key=lambda r: float(r["E_data"]))
                    t.add_row(
                        nl, ml, bl,
                        str(mv) if mv != 0 else "",
                        ibp,
                        fmt_val(float(best["E_data"])),
                        fmt_val(float(best.get("cond_theta", 0)), 1),
                    )
    return t


def build_anisotropic_summary_table(rows):
    """Equivariant vs best non-equivariant summary, anisotropic noise.

    Like ``build_standard_summary_table`` but adds a Profile column.
    """
    from .summary_metadata import LIB_NAMES, METHOD_MAP, PROFILE_ORDER, PROFILE_LABELS

    t = MarkdownTable(
        [
            "Noise", "Method", "Backend", "Profile", "Eq best",
            "Multigraphs fp/fn", "Multigraphs E",
            "Non-eq best", "Non-eq fp/fn", "Non-eq E",
        ],
        align=["left"] + ["center"] * 4 + ["left"] * 5,
    )

    for noise, noise_label in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
        for method, method_label in METHOD_MAP.items():
            for be, be_label in [("spectral", "spectral"), ("polynomial", "polynomial")]:
                for pt, pd in PROFILE_ORDER:
                    label = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")

                    # Filter rows by profile
                    eq_rows = [
                        r for r in rows
                        if r.get("profile_type", "") == pt
                        and r.get("profile_direction", "") == pd
                        and r["library"] == "multigraphs"
                        and str(r["noise"]) == str(noise)
                        and r["method"] == method
                        and r["backend"] == be
                        and r.get("fp", "nan") != "nan"
                    ]
                    ne_rows = [
                        r for r in rows
                        if r.get("profile_type", "") == pt
                        and r.get("profile_direction", "") == pd
                        and r["library"] not in ("multigraphs", "mwsindy")
                        and str(r["noise"]) == str(noise)
                        and r["method"] == method
                        and r["backend"] == be
                        and r.get("fp", "nan") != "nan"
                    ]

                    eq = best_in_group(eq_rows) if eq_rows else None
                    ne = best_in_group(ne_rows) if ne_rows else None

                    if eq is None and ne is None:
                        continue

                    winner = eq_wins(eq, ne)
                    eq_fpf = f"{eq['fp']}/{eq['fn']}{status_mark(eq)}" if eq else "\u2014"
                    eq_E = fmt_val(eq["E_coeff"]) if eq else "\u2014"

                    if ne is None:
                        ne_lib = "\u2014"
                        ne_fpf = "\u2014"
                        ne_E = "\u2014"
                    else:
                        ne_lib = "\u2014" if is_failure(ne) else (
                            f"{LIB_NAMES.get(ne['library'], '')}/"
                            f"{ALGORITHM_LABELS.get(ne.get('algorithm', ''), ne.get('algorithm', ''))}"
                        )
                        ne_fpf = f"{ne['fp']}/{ne['fn']}{status_mark(ne)}"
                        ne_E = fmt_val(ne["E_coeff"])

                    t.add_row(
                        noise_label, method_label, be_label, label, winner,
                        eq_fpf, eq_E, ne_lib, ne_fpf, ne_E,
                    )
    return t


def build_anisotropic_detailed_table(rows, method, method_label, backend, backend_label):
    """Detailed per-library table for one (method, backend) with profiles."""
    from .summary_metadata import PROFILE_ORDER, PROFILE_LABELS

    t = MarkdownTable(
        [
            "Library", "Noise", "Algorithm", "d_tol/\u03b1", "Profile",
            "fp/fn", "n_active", "E_coeff", "Time (s)", "Exact?",
        ],
    )

    for lib, lib_label in [
        ("standard", "PDE-FIND"),
        ("wsindy", "WSINDy"),
        ("polynomial_jet", "PolyJet"),
        ("multigraphs", "Multigraphs"),
    ]:
        for n, nl in [("0.0", "Clean"), ("0.02", "2%"), ("0.05", "5%")]:
            for alg, alg_label in [("strridge", "STRidge"), ("mstls", "MSTLS"), ("rfe", "RFE")]:
                for pt, pd in PROFILE_ORDER:
                    plabel = PROFILE_LABELS.get((pt, pd), f"{pt}/{pd}")
                    g = [
                        r for r in rows
                        if r.get("profile_type", "") == pt
                        and r.get("profile_direction", "") == pd
                        and r["library"] == lib
                        and str(r["noise"]) == str(n)
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
                    p = b.get("mstls_alpha", "") if alg == "mstls" else b.get("d_tol", "")
                    t.add_row(
                        lib_label, nl, alg_label, str(p), plabel,
                        fpf, str(b["n_active"]), E, tim, ex,
                    )
    return t
