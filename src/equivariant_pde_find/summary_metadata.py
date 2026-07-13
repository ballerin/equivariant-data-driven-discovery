"""
Shared constants and equation metadata for sweep summarization scripts.

All summarizers (build_summary.py, build_anisotropic_summary.py,
summarize_multigraph_comparison.py) import from this module to avoid
duplicating library labels, equation LaTeX, term lists, and the legend.
"""

# ═══════════════════════════════════════════════════════════════════════════
# Library display names
# ═══════════════════════════════════════════════════════════════════════════

LIB_NAMES = {
    "multigraphs": "Multigraphs",
    "mod_multigraph": "Modified Multigraphs",
    "standard": "PDE-FIND",
    "polynomial_jet": "PolyJet",
    "wsindy": "WSINDy",
    "mwsindy": "MWSINDy",
}

METHOD_MAP = {
    "pointwise": "pointwise",
    "weak_direct": "weak (no IBP)",
    "weak_ibp": "weak (IBP)",
}

# TODO: this can be done programmatically
NOISE_LABELS = {
    "0.0": "clean",
    "0.001": "0.1%",
    "0.005": "0.5%",
    "0.01": "1%",
    "0.02": "2%",
    "0.05": "5%",
}

ALGORITHM_LABELS = {
    "stride": "STRidge",
    "sr3": "SR3",
    "mstls": "MSTLS",
    "rfe": "RFE",
    "brute": "Brute-force",
}

BACKEND_LABELS = {
    "spectral": "Spectral derivatives",
    "polynomial": "Polynomial stencils",
}

# ═══════════════════════════════════════════════════════════════════════════
# Anisotropic noise profile labels
# ═══════════════════════════════════════════════════════════════════════════

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

# ═══════════════════════════════════════════════════════════════════════════
# Equation display names and LaTeX (for multigraph comparison summaries)
# ═══════════════════════════════════════════════════════════════════════════

EQ_NAMES = {
    "ac": "Allen-Cahn",
    "sh": "Swift-Hohenberg",
    "ch": "Cahn-Hilliard",
    "pme": "Porous Medium",
    "ks_cell": "KS (cellular)",
    "ks_wchaos": "KS (weak-chaotic)",
    "ks_chaos": "KS (chaotic)",
}

EQ_TRUTH = {
    "ac": r"$f_t = \varepsilon^2 \Delta f + f - f^3$",
    "sh": r"$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3$",
    "ch": r"$f_t = -\Delta^2 f + \Delta(f^3 - f)$",
    "pme": (
        r"$f_t = \Delta(f^3)/3$ (mod) = "
        r"$f^2\Delta f + 2f\|\nabla f\|^2$ (mg)"
    ),
    "ks_cell": r"$f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$",
    "ks_wchaos": r"$f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$",
    "ks_chaos": r"$f_t = -\Delta f - \Delta^2 f - \|\nabla f\|^2/2$",
}

EQ_TERMS = {
    "ac": ("f, f\u00b3, \u0394f", "f, f\u00b3, \u0394f"),
    "sh": (
        "f, f\u00b3, \u0394f, \u0394\u00b2f",
        "f, f\u00b3, \u0394f, \u0394\u00b2f",
    ),
    "ks_cell": (
        "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
        "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
    ),
    "ks_wchaos": (
        "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
        "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
    ),
    "ks_chaos": (
        "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
        "\u0394f, \u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
    ),
    "ch": (
        "\u0394f, f\u00b2\u00b7\u0394f, f\u00b7\u2016\u2207f\u2016\u00b2, \u0394\u00b2f",
        "\u0394f, \u0394f\u00b3, \u0394\u00b2f",
    ),
    "pme": (
        "f\u00b2\u00b7\u0394f, f\u00b7\u2016\u2207f\u2016\u00b2",
        "\u0394f\u00b3",
    ),
}

# ═══════════════════════════════════════════════════════════════════════════
# Per-equation metadata for build_summary and build_anisotropic_summary
# ═══════════════════════════════════════════════════════════════════════════

# Each equation entry: title, latex, equivariant_terms, standard_terms, notes.
# standard_terms is a list of (term, pde_find, wsindy, poly_jet, multigraph).
# equivariant_terms is a list of (term, coeff_string).

EQUATION_META = {
    "pme": {
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
            "minimum fn is 2.  "
            "WSINDy represents PME exactly via the chain rule: "
            "$f_t = \\frac13\\partial_x^2(f^3) + "
            "\\frac13\\partial_y^2(f^3)$, which "
            "expands to all 4 jet terms."
        ),
    },
    "ac": {
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
            "All four libraries can represent every term, and "
            "this is considered the simplest equation in the "
            "benchmark set."
        ),
    },
    "sh": {
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
            "giving it a sparsity advantage"
        ),
    },
    "ch": {
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
            "$\\partial_{xx}(f^2)$).  Three KS regimes are simulated: "
            "cellular ($L=4\\pi$), weak-chaotic ($L=8\\pi$), "
            "chaotic ($L=16\\pi$)."
        ),
    },
}

# ═══════════════════════════════════════════════════════════════════════════
# Legend for summary tables
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
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def get_meta(eq_name):
    """Look up EQUATION_META by equation name (prefix match)."""
    for prefix, m in EQUATION_META.items():
        if eq_name.startswith(prefix):
            return m
    return None
