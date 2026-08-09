"""
O(2)-equivariant multigraph operators on R^2, with p <= 3 (edges) and n <= 4 (nodes).

Implements all base multigraph operators N_Gamma for the required p and n as described in:

  "Equivariant nonlinear partial differential operators on
   constant curvature spaces"  - Ballerin & Grong, 2026.

TODO: Update references to the final published version of the paper once available.

Background
----------
Each operator N_Gamma[f] corresponds to a multigraph Gamma whose vertices
represent copies of f and whose edges encode metric contractions of partial
derivatives.  Every non-isolated vertex of degree d contributes a rank-d
covariant derivative tensor of f; isolated vertices contribute bare factors
of f.  Edges contract these tensors pairwise via the Euclidean metric.

Parameters
----------
p  (edges in Gamma) controls the total differential order: order = 2p.
n  (total vertices = non-isolated k + isolated m) is the polynomial degree.

The complete non-redundant basis for O(2)-equivariant polynomial nonlinear PDOs
on R^2 up to order 2p and polynomial degree n is obtained by calling
``get_operators(max_p, max_n)``.

Structure
---------
``CONNECTED_OPERATORS`` is a dict  {p: [(func, description, k, excluded_2d), ...]}
where each entry is a connected multigraph operator with p edges, k non-isolated
vertices, and excluded_2d=True for the two operators that vanish identically in R^2.

``get_operators(max_p, max_n)`` returns a flat list of ``(func, desc)`` pairs
covering all products of non-excluded connected operators with total p' <= max_p,
and all allowed polynomial multipliers f^m (m >= 0) such that k + m <= max_n.
Operators are ordered by (total_p, k, m).

2D-specific kernel elements (p=3)
----------------------------------
On R^2 two multigraph operators are linearly dependent on the others and are
excluded from the non-redundant basis.  Both vanish identically in d=2 but are
independent in d≥3.

  k=3  Triangle graph (three vertices each of degree 2, mutually connected):
       tr((∇²f)^3)  is linearly dependent via the 2-D Cayley-Hamilton identity:
         tr((∇²f)^3) = (Δf)^3 - 3(Δf)·det(∇²f),
         det(∇²f) = [(Δf)^2 - ||∇²f||^2] / 2.
       Excluded: only (Δf)^3 and ||∇²f||^2·Δf survive as k=3 basis elements.

  k=4  Path-of-4 graph (four vertices of degrees 1,2,2,1 forming a path):
       ||∇²f(∇f,·)||^2 = Σ_i (Σ_j f_{ij} f_j)^2  is linearly dependent:
         ||∇f||^2·(Δf)^2 - ||∇f||^2·||∇²f||^2
           + 2||∇²f(∇f,·)||^2 - 2Δf·∇²f(∇f,∇f) = 0   (in R^2).
       Excluded: star graph ∇³f(∇f,∇f,∇f) and the three disconnected k=4
       operators survive as basis elements.

N-dimensional generalization (max_p <= 2 only)
------------------------------------------------
Every p<=2 operator (both in this module's original basis and in the
modified/IBP-friendly basis below) is an isotropic scalar contraction of
grad/Hessian/Laplacian and therefore generalizes to an arbitrary number of
spatial dimensions.  Each is defined symbolically once, as a formula over
``jet_algebra.JetExpr`` (see that module), and evaluated numerically via
``_evaluate_jet_expr``; the same symbolic formula also drives the analytic
jet-space decomposition used by ``library/conversions.py`` for arbitrary
ndim, so the two can never silently disagree.

p=3 operators (below) are NOT generalized: three of them already require a
full third-order derivative tensor contraction and were never implemented
beyond ndim=3 in the first place.  Calling any p=3 operator with a field of
ndim not in (2, 3) raises ``NotImplementedError`` with a clear message; use
``max_p<=2`` for N-D (ndim>3) support.
"""

from functools import lru_cache

import numpy as np

from .derivatives import (
    bilaplacian,
    grad,
    hessian,
    laplacian,
    third_order,
    trilaplacian,
)
from .jet_algebra import JetExpr

# -----------------------------------------------------------------------
# Private base operator functions  (no polynomial prefactor)
# Naming: _op_<description>
# Each function has signature  f(f, Lx, Ly) -> ndarray of shape f.shape
# -----------------------------------------------------------------------

# --- p = 0 ---


def _op_const(f, Lx, Ly, Lz=None, cache=None):
    """
    p=0, k=0. Empty graph: the constant 1.

    Graph:  ∅
    """
    return np.ones_like(f)


def _op_identity(f, Lx, Ly, Lz=None, cache=None):
    """
    p=0, k=1. Single isolated vertex: f.

    Graph:  X
    """
    return f.copy()


# -----------------------------------------------------------------------
# Symbolic formulas for max_p<=2 (numeric- and jet-decomposition-generic).
# See jet_algebra.py's module docstring for the design.
# -----------------------------------------------------------------------


def _evaluate_jet_expr(expr, f, Lx, Ly, Lz, cache=None):
    """Numerically evaluate a JetExpr for field f (any ndim).

    If *cache* is provided (a dict), derivative computations are shared
    across calls.  Otherwise a local cache is created.
    """
    from . import derivatives

    Ls = derivatives.axis_lengths(f, Lx, Ly, Lz)
    if cache is None:
        cache = {}

    def get(multi_index):
        if multi_index not in cache:
            cache[multi_index] = (
                f if not any(multi_index) else derivatives.partial_nd(f, Ls, multi_index)
            )
        return cache[multi_index]

    total = np.zeros(f.shape, dtype=float)
    for mono, coeff in expr.monomials():
        term = np.full(f.shape, float(coeff))
        for factor in mono:
            term = term * get(factor)
        total = total + term
    return total


def _make_symbolic_op(sym_fn):
    """Wrap a symbolic formula ``sym_fn(ndim) -> JetExpr`` as a numeric
    operator function ``op(f, Lx, Ly, Lz=None, cache=None) -> ndarray``."""

    def op(f, Lx, Ly, Lz=None, cache=None, _sym=sym_fn):
        return _evaluate_jet_expr(_sym(f.ndim), f, Lx, Ly, Lz, cache=cache)

    return op


@lru_cache(maxsize=None)
def _sym_Delta(ndim):
    F = JetExpr.var(ndim)
    return sum((F.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_grad_sq(ndim):
    F = JetExpr.var(ndim)
    return sum((F.deriv(i) * F.deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_bilap(ndim):
    F = JetExpr.var(ndim)
    DF = sum((F.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))
    return sum((DF.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_grad_Delta_grad(ndim):
    F = JetExpr.var(ndim)
    DF = sum((F.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))
    return sum(
        (DF.deriv(i) * F.deriv(i) for i in range(ndim)), JetExpr.zero(ndim)
    )


@lru_cache(maxsize=None)
def _sym_hess_norm_sq(ndim):
    F = JetExpr.var(ndim)
    return sum(
        (F.deriv(i).deriv(j) * F.deriv(i).deriv(j)
         for i in range(ndim) for j in range(ndim)),
        JetExpr.zero(ndim),
    )


@lru_cache(maxsize=None)
def _sym_hess_grad_grad(ndim):
    F = JetExpr.var(ndim)
    return sum(
        (F.deriv(i).deriv(j) * F.deriv(i) * F.deriv(j)
         for i in range(ndim) for j in range(ndim)),
        JetExpr.zero(ndim),
    )


@lru_cache(maxsize=None)
def _sym_Delta_f2(ndim):
    """Δ(f²) = 2f·Δf + 2‖∇f‖². Replacement for f·Δf."""
    F = JetExpr.var(ndim)
    G = F * F
    return sum((G.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_Delta_f3(ndim):
    """Δ(f³) = 3f²·Δf + 6f·‖∇f‖². Replacement for f²·Δf."""
    F = JetExpr.var(ndim)
    G = F * F * F
    return sum((G.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_Delta_grad_sq(ndim):
    """Δ(‖∇f‖²) = 2⟨∇Δf,∇f⟩ + 2‖∇²f‖². Replacement for ⟨∇Δf,∇f⟩."""
    G = _sym_grad_sq(ndim)
    return sum((G.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_Delta_f_grad_sq(ndim):
    """Δ(f·‖∇f‖²). Replacement for f·⟨∇Δf,∇f⟩."""
    F = JetExpr.var(ndim)
    G = F * _sym_grad_sq(ndim)
    return sum((G.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_Delta2_f2(ndim):
    """Δ²(f²). Replacement for f·Δ²f."""
    F = JetExpr.var(ndim)
    G = F * F
    D1 = sum((G.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))
    return sum((D1.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


@lru_cache(maxsize=None)
def _sym_Delta2_f3(ndim):
    """Δ²(f³). Replacement for f²·Δ²f."""
    F = JetExpr.var(ndim)
    G = F * F * F
    D1 = sum((G.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))
    return sum((D1.deriv(i).deriv(i) for i in range(ndim)), JetExpr.zero(ndim))


# --- p = 1 ---

_op_Delta = _make_symbolic_op(_sym_Delta)
"""p=1, k=1. One vertex with a single loop: Δf = Σ_i f_ii."""

_op_grad_sq = _make_symbolic_op(_sym_grad_sq)
"""p=1, k=2. Single edge between two vertices: ‖∇f‖² = Σ_i f_i²."""


# --- p = 2 ---

_op_bilap = _make_symbolic_op(_sym_bilap)
"""p=2, k=1. One vertex with two loops: Δ²f."""

_op_grad_Delta_grad = _make_symbolic_op(_sym_grad_Delta_grad)
"""p=2, k=2. Single edge with one loop at an endpoint: ⟨∇Δf, ∇f⟩ = Σ_i (Δf)_i f_i."""

_op_hess_norm_sq = _make_symbolic_op(_sym_hess_norm_sq)
"""p=2, k=2. Double edge between two vertices: ‖∇²f‖² = Σ_ij f_ij²."""

_op_hess_grad_grad = _make_symbolic_op(_sym_hess_grad_grad)
"""p=2, k=3. Path of three vertices: ∇²f(∇f,∇f) = Σ_ij f_ij f_i f_j."""


# --- p = 2, modified multigraph operators ---
# These are the 6 operators in the IBP-friendly basis that differ from the
# original multigraph basis -- each is Δ(something) or Δ²(something),
# expanded via the product rule (the symbolic .deriv() chain rule handles
# this automatically, so no compounded numerical differentiation occurs:
# every jet monomial is still evaluated via a single joint derivative call).

_op_Delta_f2 = _make_symbolic_op(_sym_Delta_f2)
_op_Delta_f3 = _make_symbolic_op(_sym_Delta_f3)
_op_Delta_grad_sq = _make_symbolic_op(_sym_Delta_grad_sq)
_op_Delta_f_grad_sq = _make_symbolic_op(_sym_Delta_f_grad_sq)
_op_Delta2_f2 = _make_symbolic_op(_sym_Delta2_f2)
_op_Delta2_f3 = _make_symbolic_op(_sym_Delta2_f3)


# --- p = 3 ---


def _op_trilap(f, Lx, Ly, Lz=None):
    """
    p=3, k=1. One vertex with three loops: Δ³f.
    Vertex degree: 6.

    Graph:  ↺↺↺
            X
    """
    return trilaplacian(f, Lx, Ly, Lz)


def _op_hess_hess_Delta(f, Lx, Ly, Lz=None):
    """
    p=3, k=2. Double edge with one loop at an endpoint: ⟨∇²f, ∇²Δf⟩_F = Σ_{ij} f_{ij} (Δf)_{ij}.
    The looped vertex (degree 4) contributes ∇²Δf; the plain vertex (degree 2) contributes ∇²f.
    Vertex degrees: 4, 2.

    Graph:  ↺
            X===X
    """
    Df = laplacian(f, Lx, Ly, Lz)
    if f.ndim == 3:
        fxx, fxy, fxz, fyy, fyz, fzz = hessian(f, Lx, Ly, Lz)
        Dfxx, Dfxy, Dfxz, Dfyy, Dfyz, Dfzz = hessian(Df, Lx, Ly, Lz)
        return (fxx * Dfxx + fyy * Dfyy + fzz * Dfzz
                + 2 * (fxy * Dfxy + fxz * Dfxz + fyz * Dfyz))
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    Dfxx, Dfxy, Dfyy = hessian(Df, Lx, Ly)
    return fxx * Dfxx + 2 * fxy * Dfxy + fyy * Dfyy


# --- p = 3, connected operators ---


def _op_nabla3_nabla3(f, Lx, Ly, Lz=None):
    """
    p=3, k=4. Star graph: one central vertex (degree 3) with three pendant vertices (degree 1).
    ∇³f(∇f,∇f,∇f).
    The center contributes ∇³f; each leaf contributes ∇f.
    Vertex degrees: 1, 1, 1, 3 (center).

    Graph:      X
                |
            X---X---X
    """
    if f.ndim == 3:
        derivs = third_order(f, Lx, Ly, Lz)
        fx, fy, fz = grad(f, Lx, Ly, Lz)
        # In 3D the star graph gives a more complex expression.
        # For simplicity, delegate to the general 3D third-order contraction.
        raise NotImplementedError(
            "_op_nabla3_nabla3 not yet implemented for 3-D. "
            "Use max_p <= 2 for 3-D experiments."
        )
    fxxx, fxxy, fxyy, fyyy = third_order(f, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return fxxx * fx**3 + 3 * fxxy * fx**2 * fy + 3 * fxyy * fx * fy**2 + fyyy * fy**3


def _op_path_of_4(f, Lx, Ly, Lz=None):
    """
    p=3, k=4. Path of four vertices. EXCLUDED in R² via the P⁰ identity.
    ‖∇²f(∇f,·)‖² = Σ_i (Σ_j f_{ij} f_j)².
    Vertex degrees: 1, 2, 2, 1.

    Graph:  X---X---X---X
    """
    if f.ndim == 3:
        fx, fy, fz = grad(f, Lx, Ly, Lz)
        fxx, fxy, fxz, fyy, fyz, fzz = hessian(f, Lx, Ly, Lz)
        return (
            (fxx * fx + fxy * fy + fxz * fz) ** 2
            + (fxy * fx + fyy * fy + fyz * fz) ** 2
            + (fxz * fx + fyz * fy + fzz * fz) ** 2
        )
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return (fxx * fx + fxy * fy) ** 2 + (fxy * fx + fyy * fy) ** 2


def _op_triangle(f, Lx, Ly, Lz=None):
    """
    p=3, k=3. Triangle (three vertices of degree 2, each pair connected by one edge).
    EXCLUDED in R² via the Q⁰ (Cayley–Hamilton) identity.
    tr((∇²f)³).
    Vertex degrees: 2, 2, 2.

    Graph:      X
               / \
              X---X
    """
    if f.ndim == 3:
        fxx, fxy, fxz, fyy, fyz, fzz = hessian(f, Lx, Ly, Lz)
        # tr(H³) for symmetric 3×3 Hessian H = [[a,d,e],[d,b,f],[e,f,c]]:
        #   a³+b³+c³ + 3(ad²+ae²+bd²+bf²+ce²+cf²) + 6def
        return (
            fxx**3 + fyy**3 + fzz**3
            + 3.0 * (
                fxx * fxy**2 + fxx * fxz**2
                + fyy * fxy**2 + fyy * fyz**2
                + fzz * fxz**2 + fzz * fyz**2
            )
            + 6.0 * fxy * fxz * fyz
        )
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return fxx**3 + fyy**3 + 3 * fxy**2 * (fxx + fyy)


def _op_hess_inner_nabla3_grad(f, Lx, Ly, Lz=None):
    """
    p=3, k=3. Double edge between two vertices, one of which also has a single edge to a
    third vertex.
    ⟨∇²f, ∇³f(·,·,∇f)⟩_F = Σ_{ijk} f_{ij} f_{ijk} f_k.
    Vertex degrees: 2, 3, 1.

    Graph:  X===X---X
    """
    if f.ndim == 3:
        raise NotImplementedError(
            "_op_hess_inner_nabla3_grad not yet implemented for 3-D. "
            "Use max_p <= 2 for 3-D experiments."
        )
    fxxx, fxxy, fxyy, fyyy = third_order(f, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    gx = fxx * fxxx + 2 * fxy * fxxy + fyy * fxyy  # Σ_{ij} f_{ij} f_{ijx}
    gy = fxx * fxxy + 2 * fxy * fxyy + fyy * fyyy  # Σ_{ij} f_{ij} f_{ijy}
    return gx * fx + gy * fy


def _op_hess_grad_gradDelta(f, Lx, Ly, Lz=None):
    """
    p=3, k=3. Path of three vertices with one loop at an endpoint.
    ∇²f(∇f, ∇Δf) = Σ_{kl} f_{kl} f_k (∂_l Δf).
    Vertex degrees: 3, 2, 1.

    Graph:  ↺
            X---X---X
    """
    Df = laplacian(f, Lx, Ly, Lz)
    if f.ndim == 3:
        Dfx, Dfy, Dfz = grad(Df, Lx, Ly, Lz)
        fx, fy, fz = grad(f, Lx, Ly, Lz)
        fxx, fxy, fxz, fyy, fyz, fzz = hessian(f, Lx, Ly, Lz)
        return (fxx * fx * Dfx + fyy * fy * Dfy + fzz * fz * Dfz
                + fxy * (fx * Dfy + fy * Dfx)
                + fxz * (fx * Dfz + fz * Dfx)
                + fyz * (fy * Dfz + fz * Dfy))
    Dfx, Dfy = grad(Df, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return fxx * fx * Dfx + fxy * (fx * Dfy + fy * Dfx) + fyy * fy * Dfy


def _op_hessDelta_grad_grad(f, Lx, Ly, Lz=None):
    """
    p=3, k=3. Path of three vertices with one loop at the center vertex.
    ∇²Δf(∇f, ∇f) = Σ_{kl} (∂_k∂_l Δf) f_k f_l.
    Vertex degrees: 1, 4, 1.

    Graph:      ↺
            X---X---X
    """
    Df = laplacian(f, Lx, Ly, Lz)
    if f.ndim == 3:
        Dfxx, Dfxy, Dfxz, Dfyy, Dfyz, Dfzz = hessian(Df, Lx, Ly, Lz)
        fx, fy, fz = grad(f, Lx, Ly, Lz)
        return (Dfxx * fx**2 + Dfyy * fy**2 + Dfzz * fz**2
                + 2 * (Dfxy * fx * fy + Dfxz * fx * fz + Dfyz * fy * fz))
    Dfxx, Dfxy, Dfyy = hessian(Df, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return Dfxx * fx**2 + 2 * Dfxy * fx * fy + Dfyy * fy**2


def _op_nabla3_norm_sq(f, Lx, Ly, Lz=None):
    """
    p=3, k=2. Triple edge between two vertices: ‖∇³f‖².
    Both vertices contribute ∇³f (degree 3); the triple edge contracts all index triples.
    Vertex degrees: 3, 3.

    Graph:  X≡≡≡X
    """
    if f.ndim == 3:
        raise NotImplementedError(
            "_op_nabla3_norm_sq not yet implemented for 3-D. "
            "Use max_p <= 2 for 3-D experiments."
        )
    fxxx, fxxy, fxyy, fyyy = third_order(f, Lx, Ly)
    return fxxx**2 + 3 * fxxy**2 + 3 * fxyy**2 + fyyy**2


def _op_gradDelta_sq(f, Lx, Ly, Lz=None):
    """
    p=3, k=2. Single edge with one loop at each endpoint: ‖∇Δf‖².
    Each vertex (degree 3) contributes ∇Δf; the edge contracts their free indices.
    Vertex degrees: 3, 3.

    Graph:  ↺   ↺
            X---X
    """
    Df = laplacian(f, Lx, Ly, Lz)
    if f.ndim == 3:
        Dfx, Dfy, Dfz = grad(Df, Lx, Ly, Lz)
        return Dfx**2 + Dfy**2 + Dfz**2
    Dfx, Dfy = grad(Df, Lx, Ly)
    return Dfx**2 + Dfy**2


def _op_gradBilap_grad(f, Lx, Ly, Lz=None):
    """
    p=3, k=2. Single edge with two loops at one endpoint: ⟨∇Δ²f, ∇f⟩ = Σ_k (∂_k Δ²f) f_k.
    The doubly-looped vertex (degree 5) contributes ∇Δ²f; the plain vertex (degree 1)
    contributes ∇f.
    Vertex degrees: 5, 1.

    Graph:  ↺↺
            X---X
    """
    Df2 = bilaplacian(f, Lx, Ly, Lz)
    if f.ndim == 3:
        Df2x, Df2y, Df2z = grad(Df2, Lx, Ly, Lz)
        fx, fy, fz = grad(f, Lx, Ly, Lz)
        return Df2x * fx + Df2y * fy + Df2z * fz
    Df2x, Df2y = grad(Df2, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return Df2x * fx + Df2y * fy


# -----------------------------------------------------------------------
# CONNECTED_OPERATORS  -  dict  {p: [(func, desc, k, excluded_2d), ...]}
#
# k            = number of non-isolated vertices (minimum polynomial degree)
# excluded_2d  = True for operators that in R^2 are linearly dependent on others
# -----------------------------------------------------------------------

CONNECTED_OPERATORS = {
    0: [
        (_op_const, "1", 0, False),  # empty graph
        (_op_identity, "f", 1, False),  # isolated vertex
    ],
    1: [
        (_op_grad_sq, "‖∇f‖²", 2, False),  # single edge
        (_op_Delta, "Δf", 1, False),  # single loop
    ],
    2: [
        (_op_hess_grad_grad, "∇²f(∇f,∇f)", 3, False),  # path of length 2
        (_op_hess_norm_sq, "‖∇²f‖²", 2, False),  # double edge
        (_op_grad_Delta_grad, "⟨∇Δf,∇f⟩", 2, False),  # edge + loop
        (_op_bilap, "Δ²f", 1, False),  # double loop
    ],
    3: [
        (_op_path_of_4, "‖∇²f(∇f,·)‖²", 4, True),  # path-of-4 (P⁰, excluded)
        (_op_nabla3_nabla3, "∇³f(∇f,∇f,∇f)", 4, False),  # star
        (_op_triangle, "tr((∇²f)³)", 3, True),  # triangle  (Q⁰, excluded)
        (_op_hess_inner_nabla3_grad, "⟨∇²f,∇³f(·,·,∇f)⟩", 3, False),  # double-edge + single-edge
        (_op_hess_grad_gradDelta, "∇²f(∇f,∇Δf)", 3, False),  # path + side loop
        (_op_hessDelta_grad_grad, "∇²Δf(∇f,∇f)", 3, False),  # path + central loop
        (_op_nabla3_norm_sq, "‖∇³f‖²", 2, False),  # triple edge
        (_op_hess_hess_Delta, "⟨∇²f,∇²Δf⟩", 2, False),  # double edge + loop
        (_op_gradDelta_sq, "‖∇Δf‖²", 2, False),  # edge + loops at both ends
        (_op_gradBilap_grad, "⟨∇Δ²f,∇f⟩", 2, False),  # edge + double loop
        (_op_trilap, "Δ³f", 1, False),  # triple loop
    ],
}


def _require_ndim_2_or_3(func, desc):
    """Wrap a p=3 operator so it fails clearly outside ndim in (2, 3),
    instead of a confusing shape/unpack error from the underlying
    2-D/3-D-hardcoded implementation."""

    def wrapped(f, Lx, Ly, Lz=None, cache=None, _fn=func, _desc=desc):
        if f.ndim not in (2, 3):
            raise NotImplementedError(
                f"multigraph operator {_desc!r} (p=3) is only implemented for "
                f"ndim in (2, 3); got ndim={f.ndim}. Use max_p<=2 for N-D "
                f"(ndim>3) support."
            )
        return _fn(f, Lx, Ly, Lz)

    return wrapped


CONNECTED_OPERATORS[3] = [
    (_require_ndim_2_or_3(fn, desc), desc, k, excl)
    for fn, desc, k, excl in CONNECTED_OPERATORS[3]
]


# -----------------------------------------------------------------------
# Description helpers
# -----------------------------------------------------------------------

_SUPERSCRIPTS = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}


def _poly_prefix(m):
    """Return the Unicode polynomial prefix string for exponent m >= 1."""
    if m == 1:
        return "f·"
    sup = _SUPERSCRIPTS.get(m, f"^{m}")
    return f"f{sup}·"


def _make_desc(base_desc, p, m):
    """Build a human-readable description for a base operator times f^m."""
    if m == 0:
        return base_desc
    if p == 0:
        # Pure polynomial: base_desc is '1', so just show 'f', 'f²', ...
        sup = _SUPERSCRIPTS.get(m, f"^{m}")
        return "f" if m == 1 else f"f{sup}"
    return _poly_prefix(m) + base_desc


def _make_func(base_func, m):
    """Wrap base_func with an f^m polynomial factor (m >= 1)."""

    def wrapped(f, Lx, Ly, Lz=None, cache=None, _bf=base_func, _m=m):
        return f**_m * _bf(f, Lx, Ly, Lz, cache=cache)

    return wrapped


def _all_base_operators(max_p, max_n):
    """Enumerate all irredundant operators with total edge count <= max_p and
    total vertex count <= max_n.

    Generates all multisets of non-excluded connected components from
    CONNECTED_OPERATORS - including the isolated vertex (f, p=0, k=1) which
    provides the polynomial prefactors f^m - with total edge count <= max_p
    and total vertex count <= max_n.  Components are chosen in non-decreasing
    order of (p, original_index) to avoid duplicate multisets.

    Termination is guaranteed because:
    - p=0 components (isolated vertex f) do not decrease remaining_p but do
      decrease remaining_k, which is bounded below by 0.
    - p>=1 components decrease remaining_p, which is also bounded below.

    Returns
    -------
    list of (total_p, k, func, desc)
        Sorted by (total_p, k).  k is the total vertex count (non-isolated
        plus isolated), which equals the polynomial degree of the operator.
    """
    # Flat sorted list of valid components (k >= 1, non-excluded),
    # ordered by (p, original_index).  The isolated vertex 'f' (p=0, k=1)
    # comes first and accounts for polynomial prefactors f^m.
    components = []
    for p in sorted(CONNECTED_OPERATORS.keys()):
        for fn, desc, k, excl in CONNECTED_OPERATORS[p]:
            if not excl and k >= 1:
                components.append((p, k, fn, desc))
    # components[0] = (0, 1, _op_identity, 'f')  - isolated vertex
    # components[1:] = p >= 1 operators in their original order

    results = []  # (total_p, k, func, desc)

    def recurse(remaining_p, remaining_k, min_ci, n_f, base_func, base_desc, cur_k, cur_p):
        # full description: f^n_f prefix on the p>=1 part (base_desc)
        full_desc = _make_desc(base_desc, cur_p, n_f)
        # full function: f^n_f * base_func
        full_func = base_func if n_f == 0 else _make_func(base_func, n_f)
        results.append((cur_p, cur_k, full_func, full_desc))

        for ci in range(min_ci, len(components)):
            p, k, fn, desc = components[ci]
            if k > remaining_k:
                continue
            if p > 0 and p > remaining_p:
                break  # components sorted by p; all remaining p>=1 ones too expensive
            if p == 0:
                # Isolated vertex: increment f-count, p-budget unchanged
                recurse(
                    remaining_p,
                    remaining_k - k,
                    ci,
                    n_f + 1,
                    base_func,
                    base_desc,
                    cur_k + k,
                    cur_p,
                )
            else:
                # p >= 1 operator: extend the differential part
                _cf, _fn = base_func, fn
                if base_desc == "1":
                    new_base_func = fn
                    new_base_desc = desc
                else:
                    def new_base_func(f, Lx, Ly, Lz=None, cache=None, _a=_cf, _b=_fn):
                        return _a(f, Lx, Ly, Lz, cache=cache) * _b(f, Lx, Ly, Lz, cache=cache)

                    new_base_desc = base_desc + "·" + desc
                recurse(
                    remaining_p - p,
                    remaining_k - k,
                    ci,
                    n_f,
                    new_base_func,
                    new_base_desc,
                    cur_k + k,
                    cur_p + p,
                )

    recurse(max_p, max_n, 0, 0, _op_const, "1", 0, 0)
    # Canonical ordering (three-level lexicographic):
    #   1. total edge count p   - controls differential order (2p)
    #   2. polynomial degree k  - total vertex count (non-isolated + isolated)
    #   3. description string   - Unicode lexicographic order (tie-breaker)
    # This order is documented in README.md and enforced by test_multigraph_operators.py.
    results.sort(key=lambda x: (x[0], x[1], x[3]))
    return results


# -----------------------------------------------------------------------
# Public APIs
# -----------------------------------------------------------------------


def get_operators(max_p, max_n):
    """Return all equivariant operators with p <= max_p and degree <= max_n.

    Parameters
    ----------
    max_p : int  - maximum edge count (0, 1, 2, or 3); controls differential
                   order (total order = 2 * max_p).
    max_n : int  - maximum polynomial degree (= total number of f-factors,
                   counting both non-isolated vertices and extra f^m prefactors).

    Returns
    -------
    ops : list of (callable, str)
        Each entry is ``(func, description)`` where ``func(f, Lx, Ly)``
        evaluates the operator on an (Nx, Ny) field array and returns an
        array of the same shape.  Operators are ordered by (total_p, k, m).

    Examples
    --------
    >>> ops = get_operators(max_p=2, max_n=3)   # 20 operators
    >>> ops = get_operators(max_p=3, max_n=4)   # 68 operators
    """
    if max_p > max(CONNECTED_OPERATORS):
        raise ValueError(
            f"max_p={max_p} is not implemented. Available: 0..{max(CONNECTED_OPERATORS)}"
        )
    return [(func, desc) for _, _, func, desc in _all_base_operators(max_p, max_n)]


def get_operator_metadata(max_p, max_n):
    """Return ``(p, k)`` metadata for each operator in canonical order.

    Parameters
    ----------
    max_p : int
    max_n : int

    Returns
    -------
    meta : list of (int, int)
        ``(p, k)`` per operator, where ``p`` is the edge count (differential
        order = 2p) and ``k`` is the polynomial degree (total vertex count).
        Ordering matches ``get_operators(max_p, max_n)`` exactly.
    """
    if max_p > max(CONNECTED_OPERATORS):
        raise ValueError(
            f"max_p={max_p} is not implemented. Available: 0..{max(CONNECTED_OPERATORS)}"
        )
    return [(total_p, k) for total_p, k, _, _ in _all_base_operators(max_p, max_n)]


# -----------------------------------------------------------------------
# Modified Multigraph operators  (IBP-friendly basis)
# -----------------------------------------------------------------------

_MOD_MULTIGRAPH_SPEC = [
    # --- p=0 ---
    (_op_const,          "1",               0, 0),
    (_op_identity,       "f",               0, 1),
    # Pure polynomials: use direct f^m, not _make_func(_op_identity, m)
    # because _make_func multiplies the base result (f) by f^m, giving f^(m+1).
    (lambda f, Lx, Ly, Lz=None, cache=None: f**2, "f²",          0, 2),
    (lambda f, Lx, Ly, Lz=None, cache=None: f**3, "f³",          0, 3),
    # --- p=1 ---
    (_op_Delta,          "Δf",              1, 1),
    (_op_Delta_f2,       "Δf²",             1, 2),   # ← replacement
    (_op_grad_sq,        "‖∇f‖²",          1, 2),
    (_op_Delta_f3,       "Δf³",             1, 3),   # ← replacement
    (_make_func(_op_grad_sq, 1), "f·‖∇f‖²", 1, 3),
    # --- p=2 ---
    (_op_bilap,          "Δ²f",             2, 1),
    (_op_Delta2_f2,      "Δ²f²",            2, 2),   # ← replacement
    (_op_Delta_grad_sq,  "Δ‖∇f‖²",         2, 2),   # ← replacement
    (_op_hess_norm_sq,   "‖∇²f‖²",         2, 2),
    # (Δf)²: use _op_Delta square wrapped - note: original basis calls this "Δf·Δf"
    (lambda f, Lx, Ly, Lz=None, cache=None: laplacian(f, Lx, Ly, Lz)**2, "Δf·Δf", 2, 2),
    (_op_Delta2_f3,      "Δ²f³",            2, 3),   # ← replacement
    (_op_Delta_f_grad_sq,"Δ(f·‖∇f‖²)",      2, 3),   # ← replacement
    (_make_func(_op_hess_norm_sq, 1), "f·‖∇²f‖²", 2, 3),
    # f·(Δf)²: f * Δf² - original calls this "f·Δf·Δf"
    (lambda f, Lx, Ly, Lz=None, cache=None: f * laplacian(f, Lx, Ly, Lz)**2, "f·Δf·Δf", 2, 3),
    (_op_hess_grad_grad, "∇²f(∇f,∇f)",      2, 3),
    # ‖∇f‖²·Δf = ‖∇f‖² * Δf
    (lambda f, Lx, Ly, Lz=None, cache=None: (_op_grad_sq(f, Lx, Ly, Lz) * laplacian(f, Lx, Ly, Lz)),
     "‖∇f‖²·Δf", 2, 3),
]


def get_mod_multigraph_operators(max_p=2, max_n=3):
    """Return all operators in the modified multigraph (IBP-friendly) basis.

    A different basis which spans the same 20-dimensional space as the main
    multigraph basis (``get_operators``) but replaces 6 operators with
    Δ(something) or Δ²(something) forms that are better suited for
    integration by parts in the weak formulation.

    Parameters
    ----------
    max_p : int  - maximum edge count (0, 1, or 2 only)
    max_n : int  - maximum polynomial degree (0..3 only)

    Returns
    -------
    ops : list of (callable, str) operators ordered by 
            ``(total_p, k, description)``.
    """
    if max_p > 2:
        raise ValueError(
            f"max_p={max_p} not supported for modified multigraph. "
            f"Available: 0..2"
        )
    # Build {description: (p, k)} lookup for sorting.
    _pk = {desc: (p, k) for _, desc, p, k in _MOD_MULTIGRAPH_SPEC}
    results = [
        (func, desc) for func, desc, p, k in _MOD_MULTIGRAPH_SPEC
        if p <= max_p and k <= max_n
    ]
    results.sort(key=lambda x: (_pk[x[1]][0], _pk[x[1]][1], x[1]))
    return results
