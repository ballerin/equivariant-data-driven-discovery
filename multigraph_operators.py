"""
O(2)-equivariant multigraph operators on R^2, with p <= 3 (edges) and n <= 4 (nodes).

Implements all base multigraph operators N_Gamma for the required p and n as described in:

  "Equivariant nonlinear partial differential operators on
   constant curvature spaces"  — Ballerin & Grong, 2026.

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
excluded from the irredundant basis.  Both vanish identically in d=2 but are
independent in d≥3 (Theorem 2 of the paper).

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
"""

import numpy as np
from derivatives import (
    grad,
    laplacian,
    hessian,
    bilaplacian,
    trilaplacian,
    third_order,
)


# -----------------------------------------------------------------------
# Private base operator functions  (no polynomial prefactor)
# Naming: _op_<description>
# Each function has signature  f(f, Lx, Ly) -> ndarray of shape f.shape
# -----------------------------------------------------------------------

# --- p = 0 ---

def _op_const(f, Lx, Ly):
    """
    p=0, k=0. Empty graph: the constant 1.

    Graph:  ∅
    """
    return np.ones_like(f)

def _op_identity(f, Lx, Ly):
    """
    p=0, k=1. Single isolated vertex: f.

    Graph:  X
    """
    return f.copy()

# --- p = 1 ---

def _op_Delta(f, Lx, Ly):
    """
    p=1, k=1. One vertex with a single loop: Δf.
    The loop contracts both indices of ∇²f into a trace.
    Vertex degree: 2 (a loop counts as 2).

    Graph:  ↺
            X
    """
    return laplacian(f, Lx, Ly)


def _op_grad_sq(f, Lx, Ly):
    """
    p=1, k=2. Single edge between two vertices: ‖∇f‖² = f_x² + f_y².
    Each vertex contributes ∇f (degree 1); the edge contracts their free indices.
    Vertex degrees: 1, 1.

    Graph:  X---X
    """
    fx, fy = grad(f, Lx, Ly)
    return fx**2 + fy**2


# --- p = 2 ---

def _op_bilap(f, Lx, Ly):
    """
    p=2, k=1. One vertex with two loops: Δ²f.
    Vertex degree: 4.

    Graph:  ↺↺
            X
    """
    return bilaplacian(f, Lx, Ly)


def _op_grad_Delta_grad(f, Lx, Ly):
    """
    p=2, k=2. Single edge with one loop at an endpoint: ⟨∇Δf, ∇f⟩ = Σ_i (∂_i Δf)(∂_i f).
    The looped vertex (degree 3) contributes ∇Δf; the plain vertex (degree 1) contributes ∇f.
    Vertex degrees: 3, 1.

    Graph:  ↺
            X---X
    """
    Df = laplacian(f, Lx, Ly)
    Dfx, Dfy = grad(Df, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return Dfx * fx + Dfy * fy


def _op_hess_norm_sq(f, Lx, Ly):
    """
    p=2, k=2. Double edge between two vertices: ‖∇²f‖² = f_xx² + 2f_xy² + f_yy².
    Both vertices contribute ∇²f (degree 2); the double edge contracts all index pairs.
    Vertex degrees: 2, 2.

    Graph:  X===X
    """
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return fxx**2 + 2 * fxy**2 + fyy**2


def _op_hess_grad_grad(f, Lx, Ly):
    """
    p=2, k=3. Path of three vertices: ∇²f(∇f,∇f) = f_xx f_x² + 2f_xy f_x f_y + f_yy f_y².
    The center vertex (degree 2) contributes ∇²f; each endpoint (degree 1) contributes ∇f.
    Vertex degrees: 1, 2, 1.

    Graph:  X---X---X
    """
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return fxx * fx**2 + 2 * fxy * fx * fy + fyy * fy**2


# --- p = 3 ---

def _op_trilap(f, Lx, Ly):
    """
    p=3, k=1. One vertex with three loops: Δ³f.
    Vertex degree: 6.

    Graph:  ↺↺↺
            X
    """
    return trilaplacian(f, Lx, Ly)


def _op_hess_hess_Delta(f, Lx, Ly):
    """
    p=3, k=2. Double edge with one loop at an endpoint: ⟨∇²f, ∇²Δf⟩_F = Σ_{ij} f_{ij} (Δf)_{ij}.
    The looped vertex (degree 4) contributes ∇²Δf; the plain vertex (degree 2) contributes ∇²f.
    Vertex degrees: 4, 2.

    Graph:  ↺
            X===X
    """
    Df = laplacian(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    Dfxx, Dfxy, Dfyy = hessian(Df, Lx, Ly)
    return fxx * Dfxx + 2 * fxy * Dfxy + fyy * Dfyy


# --- p = 3, connected operators ---

def _op_nabla3_nabla3(f, Lx, Ly):
    """
    p=3, k=4. Star graph: one central vertex (degree 3) with three pendant vertices (degree 1).
    ∇³f(∇f,∇f,∇f) = f_xxx f_x³ + 3f_xxy f_x²f_y + 3f_xyy f_xf_y² + f_yyy f_y³.
    The center contributes ∇³f; each leaf contributes ∇f.
    Vertex degrees: 1, 1, 1, 3 (center).

    Graph:      X
                |
            X---X---X
    """
    fxxx, fxxy, fxyy, fyyy = third_order(f, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return fxxx*fx**3 + 3*fxxy*fx**2*fy + 3*fxyy*fx*fy**2 + fyyy*fy**3


def _op_path_of_4(f, Lx, Ly):
    """
    p=3, k=4. Path of four vertices. EXCLUDED in R² via the P⁰ identity:
        ‖∇f‖²·(Δf)² - ‖∇f‖²·‖∇²f‖² + 2‖∇²f(∇f,·)‖² - 2Δf·∇²f(∇f,∇f) = 0.
    ‖∇²f(∇f,·)‖² = (f_xx f_x + f_xy f_y)² + (f_xy f_x + f_yy f_y)².
    Vertex degrees: 1, 2, 2, 1.

    Graph:  X---X---X---X
    """
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return (fxx*fx + fxy*fy)**2 + (fxy*fx + fyy*fy)**2


def _op_triangle(f, Lx, Ly):
    """
    p=3, k=3. Triangle (three vertices of degree 2, each pair connected by one edge).
    EXCLUDED in R² via the Q⁰ (Cayley–Hamilton) identity:
        tr((∇²f)³) = (Δf)³ - 3Δf·det(∇²f),   det(∇²f) = [(Δf)² - ‖∇²f‖²] / 2.
    tr((∇²f)³) = f_xx³ + f_yy³ + 3f_xy²(f_xx + f_yy).
    Vertex degrees: 2, 2, 2.

    Graph:      X
               / \
              X---X
    """
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return fxx**3 + fyy**3 + 3*fxy**2*(fxx + fyy)


def _op_hess_inner_nabla3_grad(f, Lx, Ly):
    """
    p=3, k=3. Double edge between two vertices, one of which also has a single edge to a third vertex.
    ⟨∇²f, ∇³f(·,·,∇f)⟩_F = Σ_{ijk} f_{ij} f_{ijk} f_k.
    Left vertex (degree 2) contributes ∇²f; center vertex (degree 3) contributes ∇³f;
    right pendant (degree 1) contributes ∇f contracted with the free index of the center.
    Vertex degrees: 2, 3, 1.

    Graph:  X===X---X
    """
    fxxx, fxxy, fxyy, fyyy = third_order(f, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    gx = fxx*fxxx + 2*fxy*fxxy + fyy*fxyy  # Σ_{ij} f_{ij} f_{ijx}
    gy = fxx*fxxy + 2*fxy*fxyy + fyy*fyyy   # Σ_{ij} f_{ij} f_{ijy}
    return gx*fx + gy*fy


def _op_hess_grad_gradDelta(f, Lx, Ly):
    """
    p=3, k=3. Path of three vertices with one loop at an endpoint.
    ∇²f(∇f, ∇Δf) = Σ_{kl} f_{kl} f_k (∂_l Δf).
    Looped endpoint (degree 3) contributes ∇Δf; center (degree 2) contributes ∇²f;
    plain endpoint (degree 1) contributes ∇f.
    Vertex degrees: 3, 2, 1.

    Graph:  ↺
            X---X---X
    """
    Df = laplacian(f, Lx, Ly)
    Dfx, Dfy = grad(Df, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    fxx, fxy, fyy = hessian(f, Lx, Ly)
    return fxx*fx*Dfx + fxy*(fx*Dfy + fy*Dfx) + fyy*fy*Dfy


def _op_hessDelta_grad_grad(f, Lx, Ly):
    """
    p=3, k=3. Path of three vertices with one loop at the center vertex.
    ∇²Δf(∇f, ∇f) = Σ_{kl} (∂_k∂_l Δf) f_k f_l.
    Center vertex (degree 4) contributes ∇²Δf; both endpoints (degree 1) contribute ∇f.
    Vertex degrees: 1, 4, 1.

    Graph:      ↺
            X---X---X
    """
    Df = laplacian(f, Lx, Ly)
    Dfxx, Dfxy, Dfyy = hessian(Df, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return Dfxx*fx**2 + 2*Dfxy*fx*fy + Dfyy*fy**2


def _op_nabla3_norm_sq(f, Lx, Ly):
    """
    p=3, k=2. Triple edge between two vertices: ‖∇³f‖² = f_xxx² + 3f_xxy² + 3f_xyy² + f_yyy².
    Both vertices contribute ∇³f (degree 3); the triple edge contracts all index triples.
    Vertex degrees: 3, 3.

    Graph:  X≡≡≡X
    """
    fxxx, fxxy, fxyy, fyyy = third_order(f, Lx, Ly)
    return fxxx**2 + 3*fxxy**2 + 3*fxyy**2 + fyyy**2


def _op_gradDelta_sq(f, Lx, Ly):
    """
    p=3, k=2. Single edge with one loop at each endpoint: ‖∇Δf‖² = (∂_x Δf)² + (∂_y Δf)².
    Each vertex (degree 3) contributes ∇Δf; the edge contracts their free indices.
    Vertex degrees: 3, 3.

    Graph:  ↺   ↺
            X---X
    """
    Df = laplacian(f, Lx, Ly)
    Dfx, Dfy = grad(Df, Lx, Ly)
    return Dfx**2 + Dfy**2


def _op_gradBilap_grad(f, Lx, Ly):
    """
    p=3, k=2. Single edge with two loops at one endpoint: ⟨∇Δ²f, ∇f⟩ = Σ_k (∂_k Δ²f) f_k.
    The doubly-looped vertex (degree 5) contributes ∇Δ²f; the plain vertex (degree 1) contributes ∇f.
    Vertex degrees: 5, 1.

    Graph:  ↺↺
            X---X
    """
    Df2 = bilaplacian(f, Lx, Ly)
    Df2x, Df2y = grad(Df2, Lx, Ly)
    fx, fy = grad(f, Lx, Ly)
    return Df2x*fx + Df2y*fy


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# CONNECTED_OPERATORS  —  dict  {p: [(func, desc, k, excluded_2d), ...]}
#
# k            = number of non-isolated vertices (minimum polynomial degree)
# excluded_2d  = True for operators that in R^2 are linearly dependent on others
# -----------------------------------------------------------------------

CONNECTED_OPERATORS = {
    0: [
        (_op_const,                  '1',                   0, False),  # empty graph
        (_op_identity,               'f',                   1, False),  # isolated vertex
    ],
    1: [
        (_op_grad_sq,                '‖∇f‖²',               2, False),  # single edge
        (_op_Delta,                  'Δf',                  1, False),  # single loop
    ],
    2: [
        (_op_hess_grad_grad,         '∇²f(∇f,∇f)',          3, False),  # path of length 2
        (_op_hess_norm_sq,           '‖∇²f‖²',              2, False),  # double edge
        (_op_grad_Delta_grad,        '⟨∇Δf,∇f⟩',            2, False),  # edge + loop
        (_op_bilap,                  'Δ²f',                 1, False),  # double loop
    ],
    3: [
        (_op_path_of_4,              '‖∇²f(∇f,·)‖²',        4, True),   # path-of-4 (P⁰, excluded)
        (_op_nabla3_nabla3,          '∇³f(∇f,∇f,∇f)',       4, False),  # star
        (_op_triangle,               'tr((∇²f)³)',          3, True),   # triangle  (Q⁰, excluded)
        (_op_hess_inner_nabla3_grad, '⟨∇²f,∇³f(·,·,∇f)⟩',   3, False),  # double-edge + single-edge
        (_op_hess_grad_gradDelta,    '∇²f(∇f,∇Δf)',         3, False),  # path + side loop
        (_op_hessDelta_grad_grad,    '∇²Δf(∇f,∇f)',         3, False),  # path + central loop
        (_op_nabla3_norm_sq,         '‖∇³f‖²',              2, False),  # triple edge
        (_op_hess_hess_Delta,        '⟨∇²f,∇²Δf⟩',          2, False),  # double edge + loop
        (_op_gradDelta_sq,           '‖∇Δf‖²',              2, False),  # edge + loops at both ends
        (_op_gradBilap_grad,         '⟨∇Δ²f,∇f⟩',           2, False),  # edge + double loop
        (_op_trilap,                 'Δ³f',                 1, False),  # triple loop
    ],
}


# -----------------------------------------------------------------------
# Description helpers
# -----------------------------------------------------------------------

_SUPERSCRIPTS = {2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹'}


def _poly_prefix(m):
    """Return the Unicode polynomial prefix string for exponent m >= 1."""
    if m == 1:
        return 'f·'
    sup = _SUPERSCRIPTS.get(m, f'^{m}')
    return f'f{sup}·'


def _make_desc(base_desc, p, m):
    """Build a human-readable description for a base operator times f^m."""
    if m == 0:
        return base_desc
    if p == 0:
        # Pure polynomial: base_desc is '1', so just show 'f', 'f²', ...
        sup = _SUPERSCRIPTS.get(m, f'^{m}')
        return 'f' if m == 1 else f'f{sup}'
    return _poly_prefix(m) + base_desc


def _make_func(base_func, m):
    """Wrap base_func with an f^m polynomial factor (m >= 1)."""
    def wrapped(f, Lx, Ly, _bf=base_func, _m=m):
        return f**_m * _bf(f, Lx, Ly)
    return wrapped


def _all_base_operators(max_p, max_n):
    """Enumerate all irredundant operators with total edge count <= max_p and
    total vertex count <= max_n.

    Generates all multisets of non-excluded connected components from
    CONNECTED_OPERATORS — including the isolated vertex (f, p=0, k=1) which
    provides the polynomial prefactors f^m — with total edge count <= max_p
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
    # components[0] = (0, 1, _op_identity, 'f')  — isolated vertex
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
                recurse(remaining_p, remaining_k - k, ci, n_f + 1,
                        base_func, base_desc, cur_k + k, cur_p)
            else:
                # p >= 1 operator: extend the differential part
                _cf, _fn = base_func, fn
                if base_desc == '1':
                    new_base_func = fn
                    new_base_desc = desc
                else:
                    new_base_func = lambda f, Lx, Ly, _a=_cf, _b=_fn: _a(f, Lx, Ly) * _b(f, Lx, Ly)
                    new_base_desc = base_desc + '·' + desc
                recurse(remaining_p - p, remaining_k - k, ci, n_f,
                        new_base_func, new_base_desc, cur_k + k, cur_p + p)

    recurse(max_p, max_n, 0, 0, _op_const, '1', 0, 0)
    results.sort(key=lambda x: (x[0], x[1]))
    return results


# -----------------------------------------------------------------------
# Public APIs
# -----------------------------------------------------------------------

def get_operators(max_p, max_n):
    """Return all equivariant operators with p <= max_p and degree <= max_n.

    Parameters
    ----------
    max_p : int  — maximum edge count (0, 1, 2, or 3); controls differential
                   order (total order = 2 * max_p).
    max_n : int  — maximum polynomial degree (= total number of f-factors,
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