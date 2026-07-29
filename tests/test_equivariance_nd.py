"""N-dimensional equivariance tests — the paper's central claim, beyond 2-D.

``tests/test_equivariance.py`` establishes O(2)-equivariance on a square
grid using ``np.rot90``/``np.flip``.  Nothing tested that claim at ndim > 2,
even though the operator basis is now dimension-generic.  This file closes
that gap from two independent directions:

1. :class:`TestIsotropicInvariance` — a **grid-free, symbolic** check.  Each
   operator's jet decomposition is evaluated as a polynomial in randomly
   drawn symmetric tensors standing in for ``(f, ∇f, ∇²f, ∇³f, ∇⁴f)`` at a
   point, then re-evaluated after rotating every tensor by a random
   ``R ∈ SO(n)``.  An O(n)-invariant operator must give the same number.
   This runs at ndim = 2..6 and catches errors that no lattice symmetry can:
   ``Σ_{i≤j} f_ij²`` is invariant under every signed axis permutation but is
   *not* rotation-invariant, whereas the true ``‖∇²f‖² = Σ_{ij} f_ij²`` is.

2. :class:`TestHyperoctahedralEquivariance` — a **numeric, on-grid** check.
   The exact symmetry group of a periodic hypercubic grid is the
   hyperoctahedral group B_n (signed permutations of the axes), realized by
   ``np.transpose`` and an index-reversing reflection — both exact, with no
   interpolation.  Every multigraph operator must commute with these to
   machine precision at ndim = 3, 4, while coordinate-basis derivatives must
   not.  Where (1) tests the symbolic formulas, this tests the whole numeric
   path: derivative backend, multi-index bookkeeping and operator evaluation
   together.

Both are genuinely independent of ``jet_algebra``'s internals: the jet
monomials are recovered by parsing the *published* description strings from
``library/conversions.py`` with a locally written parser.
"""

import itertools
import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives  # noqa: E402
from equivariant_pde_find.library.conversions import (  # noqa: E402
    _multigraph_to_jet_decomposition,
    _mod_multigraph_to_jet_decomposition,
)
from equivariant_pde_find.multigraph_operators import (  # noqa: E402
    get_mod_multigraph_operators,
    get_operators,
)


@pytest.fixture(autouse=True)
def _spectral_backend():
    derivatives.set_backend("spectral")
    yield
    derivatives.set_backend("spectral")


# ---------------------------------------------------------------------------
# Locally written parser: jet description string -> monomial
# ---------------------------------------------------------------------------

_LABELS = ["x", "y", "z", "w", "v", "u", "t", "s"]
_SUPERSCRIPT = {"²": 2, "³": 3, "⁴": 4, "⁵": 5, "⁶": 6, "⁷": 7, "⁸": 8, "⁹": 9}


def parse_jet_description(desc, ndim):
    """``'f_x²·f_xx'`` -> ``((1,0), (1,0), (2,0))`` (a sorted monomial).

    Written from scratch rather than imported so that a naming/parsing
    divergence surfaces here instead of cancelling out.
    """
    if desc == "1":
        return ()
    factors = []
    for part in desc.split("·"):
        exponent = 1
        if part[-1] in _SUPERSCRIPT:
            exponent = _SUPERSCRIPT[part[-1]]
            part = part[:-1]
        if part == "f":
            multi_index = (0,) * ndim
        else:
            assert part.startswith("f_"), f"unparseable factor {part!r}"
            multi_index = [0] * ndim
            for ch in part[2:]:
                multi_index[_LABELS.index(ch)] += 1
            multi_index = tuple(multi_index)
        factors.extend([multi_index] * exponent)
    return tuple(sorted(factors))


# ---------------------------------------------------------------------------
# 1. Grid-free O(n)-invariance of the jet decompositions
# ---------------------------------------------------------------------------


def random_symmetric_tensor(rng, ndim, order):
    """A random fully-symmetric order-``order`` tensor on R^ndim.

    Full symmetry is the only constraint on the ``order``-th derivative
    tensor of a smooth function at a point, so these are legitimate generic
    jet values.
    """
    if order == 0:
        return rng.standard_normal(())
    t = rng.standard_normal((ndim,) * order)
    sym = np.zeros_like(t)
    for perm in itertools.permutations(range(order)):
        sym += np.transpose(t, perm)
    return sym / math.factorial(order)


def random_rotation(rng, ndim):
    """A uniformly random R in SO(ndim) (QR of a Gaussian, sign-corrected)."""
    q, r = np.linalg.qr(rng.standard_normal((ndim, ndim)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def rotate_tensor(t, R):
    """Apply R to every index of a tensor: T'_{i…} = R_{i a} … T_{a…}."""
    order = t.ndim
    if order == 0:
        return t
    out = t
    for axis in range(order):
        out = np.moveaxis(np.tensordot(R, out, axes=([1], [axis])), 0, axis)
    return out


def make_jet(rng, ndim, max_order=4):
    """``{order: symmetric tensor}`` standing in for (f, ∇f, ∇²f, ∇³f, ∇⁴f)."""
    return {k: random_symmetric_tensor(rng, ndim, k) for k in range(max_order + 1)}


def jet_value(jet, multi_index):
    """Value of the jet variable named by ``multi_index``, e.g. (2,1,0) reads
    the ∇³f component ``T[0, 0, 1]``."""
    idx = tuple(axis for axis, order in enumerate(multi_index) for _ in range(order))
    tensor = jet[len(idx)]
    return float(tensor[idx]) if idx else float(tensor)


def evaluate_decomposition(pairs, jet, ndim):
    """Evaluate a ``[(jet_desc, coeff), …]`` decomposition on a jet."""
    total = 0.0
    for desc, coeff in pairs:
        monomial = parse_jet_description(desc, ndim)
        term = coeff
        for multi_index in monomial:
            term *= jet_value(jet, multi_index)
        total += term
    return total


class TestIsotropicInvariance:
    """Every operator must be an isotropic (O(n)-invariant) function of the
    jet — this is exactly the property the equivariant basis is built to
    have, and the reason the discovered PDEs are coordinate-independent."""

    @pytest.mark.parametrize("ndim", [2, 3, 4, 5, 6])
    @pytest.mark.parametrize("which", ["multigraph", "mod_multigraph"])
    def test_all_operators_rotation_invariant(self, ndim, which):
        decomposition = (
            _multigraph_to_jet_decomposition(ndim) if which == "multigraph"
            else _mod_multigraph_to_jet_decomposition(ndim)
        )
        rng = np.random.default_rng(1234 + ndim)

        for trial in range(3):
            jet = make_jet(rng, ndim)
            R = random_rotation(rng, ndim)
            rotated = {k: rotate_tensor(t, R) for k, t in jet.items()}

            for desc, pairs in decomposition.items():
                before = evaluate_decomposition(pairs, jet, ndim)
                after = evaluate_decomposition(pairs, rotated, ndim)
                scale = max(abs(before), 1.0)
                assert abs(after - before) < 1e-9 * scale, (
                    f"{which} ndim={ndim} trial={trial}: operator {desc!r} is not "
                    f"O({ndim})-invariant ({before} -> {after})"
                )

    @pytest.mark.parametrize("ndim", [3, 4])
    def test_check_is_discriminating(self, ndim):
        """A non-isotropic expression must be *rejected* — otherwise the
        invariance check above would be vacuous."""
        rng = np.random.default_rng(7)
        jet = make_jet(rng, ndim)
        R = random_rotation(rng, ndim)
        rotated = {k: rotate_tensor(t, R) for k, t in jet.items()}

        # f_xx alone (one diagonal Hessian entry) is not rotation invariant.
        not_isotropic = [("f_xx", 1.0)]
        before = evaluate_decomposition(not_isotropic, jet, ndim)
        after = evaluate_decomposition(not_isotropic, rotated, ndim)
        assert abs(after - before) > 1e-6

        # Σ_{i≤j} f_ij²: hyperoctahedral-invariant but NOT rotation invariant.
        # (This is precisely the bug the on-grid tests below cannot detect.)
        half_hessian = []
        for i in range(ndim):
            for j in range(i, ndim):
                mi = [0] * ndim
                mi[i] += 1
                mi[j] += 1
                name = "f_" + "".join(_LABELS[a] * o for a, o in enumerate(mi) if o)
                half_hessian.append((name + "²", 1.0))
        before = evaluate_decomposition(half_hessian, jet, ndim)
        after = evaluate_decomposition(half_hessian, rotated, ndim)
        assert abs(after - before) > 1e-6


# ---------------------------------------------------------------------------
# 2. Hyperoctahedral (B_n) equivariance of the numeric operators on a grid
# ---------------------------------------------------------------------------


def reflect_axis(f, axis):
    """Exact realization of ``x_axis -> -x_axis`` on a periodic grid whose
    first sample sits at the origin: index ``j`` maps to ``(-j) mod N``."""
    return np.roll(np.flip(f, axis=axis), 1, axis=axis)


def signed_permutations(ndim, rng, count):
    """A sample of ``count`` elements of the hyperoctahedral group B_n, each
    as a ``(name, transform)`` pair.  Exact on a periodic hypercubic grid."""
    perms = list(itertools.permutations(range(ndim)))
    out = []
    seen = set()
    while len(out) < count:
        perm = perms[rng.integers(len(perms))]
        flips = tuple(int(b) for b in rng.integers(0, 2, size=ndim))
        key = (perm, flips)
        if key in seen:
            continue
        seen.add(key)

        def transform(f, _perm=perm, _flips=flips):
            g = np.transpose(f, _perm)
            for axis, flip in enumerate(_flips):
                if flip:
                    g = reflect_axis(g, axis)
            return g

        out.append((f"perm={perm} flips={flips}", transform))
    return out


class TestHyperoctahedralEquivariance:
    """O(n)-equivariance of the *numeric* operators, tested on the exact
    lattice symmetries so the comparison is machine-precision."""

    def _field(self, ndim, N, L, seed=0):
        """A smooth field with distinct, non-symmetric wavenumbers per axis,
        so axis permutations genuinely change it."""
        x = np.linspace(0.0, L, N, endpoint=False)
        grids = np.meshgrid(*([x] * ndim), indexing="ij")
        ks = [1, 2, 3, 1, 2, 3][:ndim]
        phase = sum(k * g for k, g in zip(ks, grids))
        return np.cos(phase) + 0.4 * np.sin(2.0 * grids[0] - grids[1])

    @pytest.mark.parametrize("ndim", [3, 4])
    @pytest.mark.parametrize("which", ["multigraph", "mod_multigraph"])
    def test_operators_commute_with_signed_permutations(self, ndim, which):
        N, L = 8, 2 * math.pi
        f = self._field(ndim, N, L)
        ops = get_operators(2, 3) if which == "multigraph" \
            else get_mod_multigraph_operators(2, 3)
        rng = np.random.default_rng(11 + ndim)

        for name, transform in signed_permutations(ndim, rng, 4):
            gf = transform(f)
            for func, desc in ops:
                lhs = func(gf, L, L, L if ndim == 3 else None)
                rhs = transform(func(f, L, L, L if ndim == 3 else None))
                np.testing.assert_allclose(
                    lhs, rhs, atol=1e-9, rtol=0,
                    err_msg=f"{which} ndim={ndim} operator {desc!r} broke under {name}",
                )

    @pytest.mark.parametrize("ndim", [3, 4])
    def test_field_is_not_symmetric_under_the_transforms(self, ndim):
        """Guard against a vacuous test: if the field happened to be
        invariant, equivariance would hold trivially for anything."""
        N, L = 8, 2 * math.pi
        f = self._field(ndim, N, L)
        rng = np.random.default_rng(11 + ndim)
        changed = 0
        for _, transform in signed_permutations(ndim, rng, 4):
            if not np.allclose(transform(f), f, atol=1e-12):
                changed += 1
        assert changed == 4, "test field is invariant under some transform"

    @pytest.mark.parametrize("ndim", [3, 4])
    def test_coordinate_derivative_is_not_equivariant(self, ndim):
        """The contrast case: a coordinate-basis derivative does *not*
        commute with axis permutation, which is the whole motivation for the
        equivariant basis."""
        N, L = 8, 2 * math.pi
        f = self._field(ndim, N, L)
        Ls = derivatives.axis_lengths(f, L, L, L)
        mi_x = tuple(1 if a == 0 else 0 for a in range(ndim))

        swap = list(range(ndim))
        swap[0], swap[1] = swap[1], swap[0]

        def transform(g):
            return np.transpose(g, swap)

        lhs = derivatives.partial_nd(transform(f), Ls, mi_x)
        rhs = transform(derivatives.partial_nd(f, Ls, mi_x))
        assert not np.allclose(lhs, rhs, atol=1e-6), (
            "f_x unexpectedly commuted with an axis swap"
        )
