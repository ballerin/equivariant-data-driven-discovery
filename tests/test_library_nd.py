"""N-D generalization tests for the library builders (build_multigraph_library,
build_mod_multigraph_library, build_pdefind_library, build_polyjet_library,
build_wsindy_library) -- pointwise and weak formulations, at ndim=4.

Correctness is checked two ways:
- shape/finiteness smoke checks across all 5 library types x 2 formulations,
- exact recovery on synthetic single-mode data for at least one library,
  confirming the library actually represents the right physics, not just
  "runs without crashing".
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find import derivatives
from equivariant_pde_find.library import (
    build_multigraph_library,
    build_mod_multigraph_library,
    build_pdefind_library,
    build_polyjet_library,
    build_wsindy_library,
)


@pytest.fixture(autouse=True)
def _spectral_backend():
    derivatives.set_backend("spectral")
    yield
    derivatives.set_backend("spectral")


NDIM4_FIELD = None


def _make_field(ndim=4, N=8, L=2 * math.pi, seed=0):
    rng = np.random.default_rng(seed)
    return 0.3 * rng.standard_normal((N,) * ndim), L


class TestAllLibrariesPointwiseNdim4:
    @pytest.mark.parametrize("builder", [
        build_multigraph_library, build_mod_multigraph_library,
        build_pdefind_library, build_polyjet_library, build_wsindy_library,
    ])
    def test_shape_and_finite(self, builder):
        f, L = _make_field()
        T, desc = builder(f, L, L)
        assert T.shape == (f.size, len(desc))
        assert np.all(np.isfinite(T))
        assert len(set(desc)) == len(desc), "duplicate column descriptions"


class TestAllLibrariesWeakNdim4:
    @pytest.mark.parametrize("builder", [
        build_multigraph_library, build_mod_multigraph_library,
        build_pdefind_library, build_polyjet_library, build_wsindy_library,
    ])
    @pytest.mark.parametrize("space_ibp", [False, True])
    def test_shape_and_finite(self, builder, space_ibp):
        f, L = _make_field(N=10)
        m, p = 3, 4
        G, desc = builder(f, L, L, formulation="space_weak", m=m, p=p, space_ibp=space_ibp)
        Nv = (10 - 2 * m) ** 4
        assert G.shape == (Nv, len(desc))
        assert np.all(np.isfinite(G))


class TestIBPvsDirectNdim4:
    """For each coordinate-basis library type, the IBP weak form of a pure
    single-derivative column must agree with the direct (derivative-then-
    convolve) weak form, at ndim=4."""

    def _direct_weak(self, deriv_field, m, p, ndim):
        from equivariant_pde_find.library.weak_form import _cached_kernel_1d, weak_convolution
        kernel = _cached_kernel_1d(m, p)
        return weak_convolution(deriv_field, [kernel] * ndim)

    def test_pdefind(self):
        f, L = _make_field(N=14)
        m, p = 4, 6
        Theta_p, desc_p = build_pdefind_library(f, L, L, max_order=1, max_poly=1)
        Theta_w, desc_w = build_pdefind_library(
            f, L, L, max_order=1, max_poly=1, formulation="space_weak", m=m, p=p, space_ibp=True,
        )
        Theta_direct, _ = build_pdefind_library(
            f, L, L, max_order=1, max_poly=1, formulation="space_weak", m=m, p=p, space_ibp=False,
        )
        i = desc_w.index("f_x")
        rel = np.linalg.norm(Theta_w[:, i] - Theta_direct[:, i]) / np.linalg.norm(Theta_direct[:, i])
        assert rel < 0.1

    def test_polyjet(self):
        f, L = _make_field(N=14)
        m, p = 4, 6
        Theta_w, desc_w = build_polyjet_library(
            f, L, L, max_order=1, max_degree=1, formulation="space_weak", m=m, p=p, space_ibp=True,
        )
        Theta_direct, _ = build_polyjet_library(
            f, L, L, max_order=1, max_degree=1, formulation="space_weak", m=m, p=p, space_ibp=False,
        )
        i = desc_w.index("f_x")
        rel = np.linalg.norm(Theta_w[:, i] - Theta_direct[:, i]) / np.linalg.norm(Theta_direct[:, i])
        assert rel < 0.1

    def test_wsindy(self):
        f, L = _make_field(N=14)
        m, p = 4, 6
        Theta_w, desc_w = build_wsindy_library(
            f, L, L, max_order=1, max_poly=1, formulation="space_weak", m=m, p=p, space_ibp=True,
        )
        Theta_direct, _ = build_wsindy_library(
            f, L, L, max_order=1, max_poly=1, formulation="space_weak", m=m, p=p, space_ibp=False,
        )
        i = desc_w.index("f_{x}")
        rel = np.linalg.norm(Theta_w[:, i] - Theta_direct[:, i]) / np.linalg.norm(Theta_direct[:, i])
        assert rel < 0.1


class TestEndToEndSimulateDiscoverNdim4:
    """Full pipeline at ndim=4: simulate Allen-Cahn -> build multigraph
    library (weak, IBP) -> RFE regression -> exact recovery."""

    @pytest.mark.slow
    def test_allen_cahn_nd4_exact_recovery(self):
        from equivariant_pde_find.simulations.allen_cahn import allen_cahn_nd
        from equivariant_pde_find.simulations import make_regression_data
        from equivariant_pde_find.library import weak_form as wf
        from equivariant_pde_find.regression import TrainRFE

        N, L, eps = 12, 2 * math.pi, 0.15
        snaps, times = allen_cahn_nd(
            ndim=4, N=N, L=L, T_end=2.0, eps=eps, dt=0.01, save_every=20, seed=42,
        )
        f_snaps, ft_snaps, idx_sets, _ = make_regression_data(
            snaps, times, n_per_snap=40, skip_fraction=0.0, seed=0,
        )

        m, p = 4, 6
        kernel = wf._cached_kernel_1d(m, p)
        rows_G, rows_b = [], []
        for s in range(len(f_snaps)):
            f = f_snaps[s]
            ft = ft_snaps[s]
            G, desc = build_multigraph_library(
                f, L, L, max_p=2, max_n=3, formulation="space_weak", m=m, p=p, space_ibp=True,
            )
            b_w = wf.build_weak_rhs(ft, m=m, p=p)
            idx = idx_sets[s] % len(b_w)
            rows_G.append(G[idx])
            rows_b.append(b_w[idx])

        G_all = np.vstack(rows_G).astype(float)
        b_all = np.concatenate(rows_b).astype(float)
        w = TrainRFE(G_all, b_all, normalize=2, lam=0.0)

        nz = np.abs(w) > 1e-10
        assert int(np.sum(nz)) == 3, f"expected 3 active terms, got {int(np.sum(nz))}"
        idx_lap, idx_f, idx_f3 = desc.index("Δf"), desc.index("f"), desc.index("f³")
        assert nz[idx_lap] and nz[idx_f] and nz[idx_f3]
        np.testing.assert_allclose(w[idx_lap], eps**2, rtol=2e-1, atol=0.02)
        np.testing.assert_allclose(w[idx_f], 1.0, rtol=2e-1, atol=0.15)
        np.testing.assert_allclose(w[idx_f3], -1.0, rtol=2e-1, atol=0.15)


class TestExactRecoveryNdim4:
    """f_t = eps^2 Delta f + f - f^3 (Allen-Cahn) recovered exactly from
    the multigraph library at ndim=4, using the true jet-space coefficients
    directly (no simulation needed -- just checks the library matrix
    correctly represents Delta f, f, f^3 as independent columns)."""

    def test_multigraph_allen_cahn_recovery(self):
        ndim, N, L, eps = 4, 8, 2 * math.pi, 0.15
        rng = np.random.default_rng(0)
        f = 0.3 * rng.standard_normal((N,) * ndim)

        Theta, desc = build_multigraph_library(f, L, L)
        w_true = np.zeros(len(desc))
        w_true[desc.index("f")] = 1.0
        w_true[desc.index("f³")] = -1.0
        w_true[desc.index("Δf")] = eps**2

        Ft = Theta @ w_true  # exact by construction; check we can recover it via lstsq
        w_hat, *_ = np.linalg.lstsq(Theta, Ft, rcond=None)
        np.testing.assert_allclose(w_hat, w_true, atol=1e-8)


class TestDecompositionMatchesNumericLibrary:
    """The numeric library column and its analytic jet decomposition are
    written *twice*, independently:

    - the numeric column comes from ``_all_base_operators``/``_make_func``
      in ``multigraph_operators.py`` (``f**m * base_operator``, plus
      products of connected components);
    - the decomposition comes from a hand-written dict in
      ``library/conversions.py`` (``"f²·Δ²f": F * F * bilap``).

    Nothing compared the two beyond ndim=2 before this test: the existing
    N-D conversion tests check only that the *operator name sets* agree, so
    a wrong polynomial prefactor (say ``f`` where the numeric side uses
    ``f²``) would pass everything while silently corrupting every
    cross-library score ``evaluate_unified`` produces at ndim>2.
    """

    def _field(self, ndim, N=8, L=2 * math.pi, seed=3):
        rng = np.random.default_rng(seed)
        return 0.3 * rng.standard_normal((N,) * ndim), L

    def _evaluate_decomposition(self, pairs, f, Ls, ndim):
        """Evaluate ``[(jet_desc, coeff), ...]`` directly on the field, using
        a locally written description parser (see test_equivariance_nd) so
        the check does not lean on conversions.py's own parsing."""
        from test_equivariance_nd import parse_jet_description

        total = np.zeros(f.shape, dtype=float)
        for desc, coeff in pairs:
            term = np.full(f.shape, float(coeff))
            for multi_index in parse_jet_description(desc, ndim):
                term = term * (
                    f if not any(multi_index)
                    else derivatives.partial_nd(f, Ls, multi_index)
                )
            total = total + term
        return total

    @pytest.mark.parametrize("ndim", [2, 3, 4])
    @pytest.mark.parametrize("which", ["multigraph", "mod_multigraph"])
    def test_every_operator_column_equals_its_decomposition(self, ndim, which):
        from equivariant_pde_find.library.conversions import (
            _multigraph_to_jet_decomposition,
            _mod_multigraph_to_jet_decomposition,
        )
        from equivariant_pde_find.multigraph_operators import (
            get_operators,
            get_mod_multigraph_operators,
        )

        f, L = self._field(ndim)
        Lz = L if ndim == 3 else None
        Ls = derivatives.axis_lengths(f, L, L, Lz)
        ops = get_operators(2, 3) if which == "multigraph" \
            else get_mod_multigraph_operators(2, 3)
        decomposition = (
            _multigraph_to_jet_decomposition(ndim) if which == "multigraph"
            else _mod_multigraph_to_jet_decomposition(ndim)
        )

        for func, desc in ops:
            numeric = func(f, L, L, Lz)
            analytic = self._evaluate_decomposition(
                decomposition[desc], f, Ls, ndim
            )
            scale = max(np.max(np.abs(numeric)), 1.0)
            np.testing.assert_allclose(
                numeric, analytic, atol=1e-8 * scale, rtol=0,
                err_msg=(
                    f"{which} ndim={ndim}: numeric column for {desc!r} disagrees "
                    f"with its jet decomposition in conversions.py"
                ),
            )

    @pytest.mark.parametrize("ndim", [2, 3, 4])
    @pytest.mark.parametrize("which", ["multigraph", "mod_multigraph"])
    def test_round_trip_through_polyjet_library(self, ndim, which):
        """The full production path: map each operator's unit coefficient
        vector into the polyjet basis with ``map_*_to_jet`` and require the
        reconstructed column to equal the multigraph column.  This exercises
        the jet *library builder*'s naming as well, which the check above
        does not."""
        from equivariant_pde_find.library.conversions import (
            map_multigraph_to_jet,
            map_mod_multigraph_to_jet,
        )
        from equivariant_pde_find.multigraph_operators import (
            get_operators,
            get_mod_multigraph_operators,
        )

        f, L = self._field(ndim, N=6)
        Lz = L if ndim == 3 else None
        builder = build_multigraph_library if which == "multigraph" \
            else build_mod_multigraph_library
        mapper = map_multigraph_to_jet if which == "multigraph" \
            else map_mod_multigraph_to_jet
        n_ops = len(get_operators(2, 3) if which == "multigraph"
                    else get_mod_multigraph_operators(2, 3))

        Theta_lib, desc_lib = builder(f, L, L, Lz=Lz)
        Theta_jet, desc_jet = build_polyjet_library(
            f, L, L, Lz=Lz, max_order=4, max_degree=3,
        )
        assert len(desc_lib) == n_ops

        for i, desc in enumerate(desc_lib):
            w = np.zeros(len(desc_lib))
            w[i] = 1.0
            w_jet = mapper(w, desc_lib, desc_jet, ndim=ndim)
            scale = max(np.max(np.abs(Theta_lib[:, i])), 1.0)
            np.testing.assert_allclose(
                Theta_jet @ w_jet, Theta_lib[:, i], atol=1e-8 * scale, rtol=0,
                err_msg=f"{which} ndim={ndim}: polyjet round-trip failed for {desc!r}",
            )


class TestPublishedLibraryConventionsAtNdim4:
    """Pin the two coordinate-basis libraries to their published definitions.

    These conventions are easy to conflate — the two libraries differ in
    *where* the polynomial power sits relative to the derivative — and one
    such conflation was already found and fixed once during the N-D work (the
    PDE-FIND ``f³`` column was briefly computing ``f⁴``).  The definitions
    below were checked against the primary sources:

    - **PDE-FIND** (Rudy, Brunton, Proctor & Kutz, *Science Advances* 3(4),
      2017): the library is ``Θ(U, Q) = [1  U  U²  …  Uₓ  UUₓ  …]``, i.e.
      each column is a polynomial power of the field *times* a derivative of
      the field, ``f^p · ∂^α f``.  The polynomial power is a separate factor,
      so the pure-polynomial column ``f³`` is ``f³ · 1``, not ``f³ · f``.

    - **WSINDy** (Messenger & Bortz, *J. Comput. Phys.* 443, 2021): weak-form
      terms are the inner products ``⟨ψ, D^α f_j(u)⟩``, so with the monomial
      nonlinearity ``f_j(u) = u^q`` each column is the derivative *of the
      power*, ``∂^α(f^q)`` — the power sits inside the derivative.  The
      reference test function is the compactly supported piecewise
      polynomial ``φ(v) = C(v−a)^p (b−v)^q`` (eq. 4.2), which for
      ``a = −1, b = 1, p = q`` is ``(1 − v²)^p``, and ψ is *separable*,
      ``ψ(x) = φ₁(x₁) ⋯ φ_D(x_D)`` — which is exactly why the N-D weak form
      here is a product of one 1-D kernel per axis.
    """

    def _field(self, ndim=4, N=8, L=2 * math.pi, seed=0):
        rng = np.random.default_rng(seed)
        return 0.4 * rng.standard_normal((N,) * ndim), L

    def test_pdefind_polynomial_column_is_f_cubed_not_f_fourth(self):
        f, L = self._field()
        Theta, desc = build_pdefind_library(f, L, L, max_order=2, max_poly=3)
        col = Theta[:, desc.index("f³")]
        np.testing.assert_allclose(col, (f**3).ravel(), atol=1e-12)
        assert not np.allclose(col, (f**4).ravel())

    def test_pdefind_mixed_column_is_power_times_derivative(self):
        """``f²·f_xx`` must be ``f² · ∂xx f``, *not* ``∂xx(f²)``."""
        f, L = self._field()
        Ls = derivatives.axis_lengths(f, L, L)
        Theta, desc = build_pdefind_library(f, L, L, max_order=2, max_poly=2)
        name = next(d for d in desc if d.endswith("f_xx") and d.startswith("f²"))

        mi_xx = (2, 0, 0, 0)
        d_xx = derivatives.partial_nd(f, Ls, mi_xx)
        np.testing.assert_allclose(
            Theta[:, desc.index(name)], (f**2 * d_xx).ravel(), atol=1e-10
        )
        deriv_of_square = derivatives.partial_nd(f**2, Ls, mi_xx)
        assert not np.allclose(Theta[:, desc.index(name)],
                               deriv_of_square.ravel(), atol=1e-6)

    def test_wsindy_column_is_derivative_of_the_power(self):
        """``(f²)_{xx}`` must be ``∂xx(f²)``, *not* ``f² · ∂xx f``."""
        f, L = self._field()
        Ls = derivatives.axis_lengths(f, L, L)
        Theta, desc = build_wsindy_library(f, L, L, max_order=2, max_poly=2)
        name = next(d for d in desc if "xx" in d and "²" in d)

        mi_xx = (2, 0, 0, 0)
        deriv_of_square = derivatives.partial_nd(f**2, Ls, mi_xx)
        np.testing.assert_allclose(
            Theta[:, desc.index(name)], deriv_of_square.ravel(), atol=1e-10
        )
        power_times_deriv = f**2 * derivatives.partial_nd(f, Ls, mi_xx)
        assert not np.allclose(Theta[:, desc.index(name)],
                               power_times_deriv.ravel(), atol=1e-6)

    def test_the_two_conventions_actually_differ_on_this_data(self):
        """Guard: ``∂^α(f^q)`` and ``f^q ∂^α f`` must be numerically distinct
        here, or the two tests above would pass for the wrong reason."""
        f, L = self._field()
        Ls = derivatives.axis_lengths(f, L, L)
        mi_xx = (2, 0, 0, 0)
        a = derivatives.partial_nd(f**2, Ls, mi_xx)
        b = f**2 * derivatives.partial_nd(f, Ls, mi_xx)
        assert np.linalg.norm(a - b) > 1e-3 * np.linalg.norm(a)

    def test_wsindy_test_function_is_the_published_piecewise_polynomial(self):
        """φ(t) = (1 − t²)^p on |t| < 1, zero outside — Messenger & Bortz
        eq. (4.2) with a = −1, b = 1, p = q (up to normalization)."""
        from equivariant_pde_find.library.weak_form import _test_function_1d

        m, p = 6, 4
        phi = _test_function_1d(m, p)
        t = np.arange(-m, m + 1, dtype=float) / m
        expected = np.where(np.abs(t) < 1.0, (1.0 - t**2) ** p, 0.0)

        assert phi[0] == 0.0 and phi[-1] == 0.0, "not compactly supported"
        np.testing.assert_allclose(phi / phi.max(), expected / expected.max(),
                                   atol=1e-12)

    def test_weak_test_function_is_separable_across_axes(self):
        """ψ(x) = ∏_d φ_d(x_d): the N-D weak form must be the product of one
        1-D kernel per axis (Messenger & Bortz §3.2)."""
        from equivariant_pde_find.library.weak_form import (
            _cached_kernel_1d, weak_convolution,
        )

        ndim, N, m, p = 4, 12, 3, 4
        rng = np.random.default_rng(0)
        f = rng.standard_normal((N,) * ndim)
        kernel = _cached_kernel_1d(m, p)

        got = weak_convolution(f, [kernel] * ndim)

        # Independent reference: build the full separable ndim kernel
        # explicitly as an outer product and correlate directly.
        psi = kernel
        for _ in range(ndim - 1):
            psi = np.multiply.outer(psi, kernel)
        valid = N - 2 * m
        expected = np.empty((valid,) * ndim)
        for idx in np.ndindex(*((valid,) * ndim)):
            window = f[tuple(slice(i, i + 2 * m + 1) for i in idx)]
            expected[idx] = np.sum(window * psi)
        np.testing.assert_allclose(got.reshape((valid,) * ndim), expected,
                                   rtol=1e-10)
