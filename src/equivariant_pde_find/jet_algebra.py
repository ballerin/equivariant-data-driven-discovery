"""
Small symbolic computer-algebra layer: polynomials in the partial
derivatives of a single scalar field f, over an ndim-dimensional domain.

Purpose
-------
Each multigraph/modified-multigraph operator (``multigraph_operators.py``)
is an isotropic tensor contraction of grad/Hessian/Laplacian, e.g.

    Delta(F, ndim) = sum(F.deriv(i).deriv(i) for i in range(ndim))
    grad_sq(F, ndim) = sum(F.deriv(i) * F.deriv(i) for i in range(ndim))

Writing operators this way (as ``.deriv(axis)``/``+``/``*`` formulas) lets
the *same* Python code run in two modes, depending on what ``F`` is:

  - **Symbolic** (``F = JetExpr.var(ndim)``): produces the exact analytic
    jet-space decomposition of the operator for arbitrary ndim (what
    ``library/conversions.py`` needs for cross-library scoring) — this
    replaces the old hand-written, 2-D-only decomposition dicts.
  - **Numeric** (``F`` wraps a real ndarray field via a derivative
    backend): produces the actual Theta/G library column.

Because both come from the *same* formula, the numeric operator and its
jet decomposition can never silently disagree.

Representation
---------------
A ``JetExpr`` is ``{monomial: coeff}``. A monomial is a sorted tuple of
"factors"; each factor is itself a multi-index (tuple of ``ndim`` nonneg
ints) representing one partial derivative of f. E.g. in 3-D, the monomial
for ``f_x * f_zz`` is ``((1, 0, 0), (0, 0, 2))``, sorted so equal monomials
always compare equal regardless of construction order. The all-zero
multi-index ``(0,)*ndim`` represents ``f`` itself; the empty monomial
``()`` represents the constant ``1``.
"""

from collections import defaultdict

# Axis labels matching equations.py::_AXIS_LABELS, extended to ndim<=8.
AXIS_LABELS = ["x", "y", "z", "w", "v", "u", "t", "s"]

_SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵",
        6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _sup(n):
    return _SUP.get(n, f"^{n}")


class JetExpr:
    """Symbolic polynomial in the partial derivatives of a scalar field f."""

    __slots__ = ("ndim", "terms")

    def __init__(self, ndim, terms=None):
        self.ndim = ndim
        self.terms = defaultdict(float)
        if terms:
            for mono, coeff in terms.items():
                self.terms[tuple(sorted(mono))] += coeff

    @classmethod
    def zero(cls, ndim):
        return cls(ndim)

    @classmethod
    def constant(cls, ndim, c):
        e = cls(ndim)
        if c != 0:
            e.terms[()] = float(c)
        return e

    @classmethod
    def var(cls, ndim, multi_index=None):
        """The base field f (multi_index=None -> all zeros), or a specific
        jet variable if multi_index is given."""
        if multi_index is None:
            multi_index = (0,) * ndim
        multi_index = tuple(multi_index)
        if len(multi_index) != ndim:
            raise ValueError(f"multi_index length {len(multi_index)} != ndim {ndim}")
        e = cls(ndim)
        e.terms[(multi_index,)] = 1.0
        return e

    def deriv(self, axis):
        """d/dx_axis of this expression (product rule + chain rule)."""
        if not (0 <= axis < self.ndim):
            raise ValueError(f"axis {axis} out of range for ndim {self.ndim}")
        result = JetExpr(self.ndim)
        for mono, coeff in self.terms.items():
            for i, factor in enumerate(mono):
                new_factor = list(factor)
                new_factor[axis] += 1
                new_mono = mono[:i] + (tuple(new_factor),) + mono[i + 1:]
                result.terms[tuple(sorted(new_mono))] += coeff
        return result

    def _coerce(self, other):
        if isinstance(other, JetExpr):
            return other
        return JetExpr.constant(self.ndim, other)

    def __add__(self, other):
        other = self._coerce(other)
        result = JetExpr(self.ndim, dict(self.terms))
        for mono, coeff in other.terms.items():
            result.terms[mono] += coeff
        return result

    __radd__ = __add__

    def __neg__(self):
        return JetExpr(self.ndim, {m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return JetExpr(self.ndim, {m: c * other for m, c in self.terms.items()})
        other = self._coerce(other)
        result = JetExpr(self.ndim)
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                result.terms[tuple(sorted(m1 + m2))] += c1 * c2
        return result

    __rmul__ = __mul__

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError("JetExpr power must be a nonnegative int")
        result = JetExpr.constant(self.ndim, 1.0)
        base = self
        for _ in range(n):
            result = result * base
        return result

    def clean(self, tol=1e-9):
        """Drop near-zero coefficients (in place) and return self."""
        self.terms = defaultdict(
            float, {m: c for m, c in self.terms.items() if abs(c) > tol}
        )
        return self

    def monomials(self):
        """Return sorted list of (monomial, coeff) for nonzero terms."""
        self.clean()
        return sorted(self.terms.items(), key=lambda mc: _monomial_sort_key(mc[0]))


def _monomial_sort_key(monomial):
    """Sort key matching the historical description convention: by total
    differential order, then lexicographically."""
    return (sum(sum(f) for f in monomial), monomial)


def factor_name(multi_index):
    """Human-readable name for a single jet variable, e.g. (2,0)->'f_xx',
    (0,)*ndim -> 'f'."""
    if not any(multi_index):
        return "f"
    parts = []
    for axis, a in enumerate(multi_index):
        if a:
            if axis >= len(AXIS_LABELS):
                raise ValueError(f"ndim too large for axis labels (axis={axis})")
            parts.append(AXIS_LABELS[axis] * a)
    return "f_" + "".join(parts)


def describe_monomial(monomial):
    """Human-readable description of a monomial, matching the existing
    convention (e.g. ((1,0),(1,0),(2,0)) -> 'f_x²·f_xx')."""
    if not monomial:
        return "1"
    counts = {}
    for factor in monomial:
        counts[factor] = counts.get(factor, 0) + 1
    parts = []
    for factor in sorted(counts, key=lambda mi: (sum(mi), mi)):
        exp = counts[factor]
        name = factor_name(factor)
        parts.append(name if exp == 1 else f"{name}{_sup(exp)}")
    return "·".join(parts)


def expr_to_dict(expr):
    """Convert a JetExpr to {description: coeff} for nonzero terms."""
    out = {}
    for mono, coeff in expr.monomials():
        out[describe_monomial(mono)] = coeff
    return out


# ---------------------------------------------------------------------------
# Reverse parser: jet monomial description string → exponent dict
# ---------------------------------------------------------------------------

_SUP_REV = {"²": 2, "³": 3, "⁴": 4, "⁵": 5, "⁶": 6, "⁷": 7, "⁸": 8, "⁹": 9}


def parse_jet_monomial(desc, ndim):
    """Parse a jet monomial description into an exponent dict.

    Inverse of ``_describe_exponents`` (polyjet.py) / ``describe_monomial``
    (this module)::

        parse_jet_monomial("f³·f_x", 2) → {(0,0): 3, (1,0): 1}
        parse_jet_monomial("1", 2)       → {}
        parse_jet_monomial("f", 2)       → {(0,0): 1}

    Parameters
    ----------
    desc : str
        A jet monomial description in the canonical format produced by
        ``describe_monomial`` or ``_describe_exponents``.  Examples:
        ``"1"``, ``"f"``, ``"f_x"``, ``"f_x²"``, ``"f³·f_x"``,
        ``"f·f_xx"``, ``"f_x²·f_xx"``.
    ndim : int
        Number of spatial dimensions (must be ≤ len(AXIS_LABELS)).

    Returns
    -------
    exponents : dict
        ``{multi_index: exponent}`` where each multi_index is an ndim-tuple
        of nonnegative ints and exponent is a positive int.  Empty dict
        for the constant ``"1"``.

    Raises
    ------
    ValueError
        If *desc* cannot be parsed.
    """
    if ndim > len(AXIS_LABELS):
        raise ValueError(f"ndim {ndim} exceeds available axis labels ({len(AXIS_LABELS)})")
    if desc == "1":
        return {}
    if not desc:
        raise ValueError(f"Empty description string")

    exponents = {}
    for part in desc.split("·"):
        part = part.strip()
        if not part:
            continue
        # Extract superscript exponent, if any
        exp = 1
        if part[-1] in _SUP_REV:
            exp = _SUP_REV[part[-1]]
            part = part[:-1]
        # Parse the factor name
        if part == "f":
            mi = (0,) * ndim
        elif part.startswith("f_"):
            mi = [0] * ndim
            for ch in part[2:]:
                try:
                    axis = AXIS_LABELS.index(ch)
                except ValueError:
                    raise ValueError(
                        f"Unknown axis label {ch!r} in {desc!r} "
                        f"(ndim={ndim}, labels={AXIS_LABELS[:ndim]})"
                    )
                if axis >= ndim:
                    raise ValueError(
                        f"Axis {ch!r} is out of range for ndim={ndim}"
                    )
                mi[axis] += 1
            mi = tuple(mi)
        else:
            raise ValueError(
                f"Cannot parse factor {part!r} in jet monomial {desc!r}. "
                f"Expected 'f' or 'f_<axes>' (e.g. 'f_x', 'f_xxy')."
            )
        exponents[mi] = exponents.get(mi, 0) + exp

    return exponents
