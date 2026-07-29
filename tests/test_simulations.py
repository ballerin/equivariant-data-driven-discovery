"""Tests for simulations.py PDE simulators — smoke tests and physical correctness.

Covers:
- All 5 PDE simulators (Allen-Cahn, Swift-Hohenberg, KS, Cahn-Hilliard, PME)
- make_regression_data (shapes, skip fraction, finite-difference accuracy)
- Physical correctness: energy dissipation, mass conservation, non-negativity
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations import (
    allen_cahn,
    cahn_hilliard,
    kuramoto_sivashinsky,
    make_regression_data,
    porous_medium,
    swift_hohenberg,
)


# ======================================================================
# make_regression_data
# ======================================================================


class TestMakeRegressionData:
    def test_output_shapes(self):
        N, Nt = 32, 10
        snaps = np.random.default_rng(0).standard_normal((Nt, N, N))
        times = np.arange(Nt) * 0.1
        f_snaps, ft_snaps, idx_sets, dt = make_regression_data(
            snaps, times, n_per_snap=100, skip_fraction=0.0, seed=0,
        )
        assert f_snaps.shape == (Nt - 2, N, N)
        assert ft_snaps.shape == (Nt - 2, N, N)
        assert len(idx_sets) == Nt - 2
        assert all(len(idx) == 100 for idx in idx_sets)
        assert dt == pytest.approx(0.1)

    def test_skip_fraction(self):
        N, Nt = 32, 20
        snaps = np.random.default_rng(0).standard_normal((Nt, N, N))
        times = np.arange(Nt) * 0.1
        f_snaps, ft_snaps, idx_sets, dt = make_regression_data(
            snaps, times, n_per_snap=50, skip_fraction=0.5, seed=0,
        )
        # skip_fraction=0.5 → keep last 10 snapshots → 8 interior
        assert f_snaps.shape[0] <= 8

    def test_finite_difference_accuracy(self):
        """Centered finite difference should be exact for a quadratic in t."""
        N, Nt = 32, 5
        dt = 0.1
        times = np.arange(Nt) * dt
        # Linear in time: f(t) = t * g(x,y) → f_t = g(x,y)
        g = np.random.default_rng(42).standard_normal((N, N))
        snaps = np.array([t * g for t in times])
        f_snaps, ft_snaps, _, _ = make_regression_data(
            snaps, times, n_per_snap=100, skip_fraction=0.0, seed=0,
        )
        np.testing.assert_allclose(ft_snaps[0], g, rtol=1e-12, atol=1e-12)

    def test_reproducibility(self):
        """Same seed → same output."""
        N, Nt = 32, 10
        snaps = np.random.default_rng(0).standard_normal((Nt, N, N))
        times = np.arange(Nt) * 0.1
        f1, ft1, _, _ = make_regression_data(
            snaps, times, n_per_snap=50, skip_fraction=0.0, seed=42,
        )
        f2, ft2, _, _ = make_regression_data(
            snaps, times, n_per_snap=50, skip_fraction=0.0, seed=42,
        )
        np.testing.assert_allclose(f1, f2)
        np.testing.assert_allclose(ft1, ft2)


# ======================================================================
# Allen-Cahn: f_t = ε²·Δf + f − f³
# ======================================================================


class TestAllenCahn:
    def test_runs_and_produces_output(self):
        N, L = 16, 2 * np.pi
        snaps, times = allen_cahn(
            N=N, L=L, T_end=0.01, eps=0.3, dt=0.005, save_every=2, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert snaps.shape[2] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))

    def test_energy_decays(self):
        """Ginzburg-Landau energy ∫[ε²|∇f|²/2 + (f²−1)²/4] dx should
        decrease (for a short simulation starting from small amplitude)."""
        N, L = 32, 2 * np.pi
        eps = 0.15
        snaps, times = allen_cahn(
            N=N, L=L, T_end=1.0, eps=eps, dt=0.01, save_every=10, seed=42,
        )
        # Compute Ginzburg-Landau energy at first and last snapshot
        dx = L / N
        def energy(f):
            from equivariant_pde_find.derivatives import spectral as sd
            fx, fy = sd.spectral_grad(f, L, L)
            grad_sq = fx ** 2 + fy ** 2
            potential = (f ** 2 - 1) ** 2 / 4
            return np.sum(eps ** 2 * grad_sq / 2 + potential) * dx ** 2

        E0 = energy(snaps[0])
        E_end = energy(snaps[-1])
        # Energy should decrease (or stay approximately the same for very
        # short simulations).  Allow a small tolerance for numerical drift.
        assert E_end <= E0 * 1.01, (
            f"Energy increased: {E0:.4f} → {E_end:.4f}"
        )

    def test_bounded_solution(self):
        """Allen-Cahn solutions should remain in [-1, 1] approximately."""
        N, L = 32, 2 * np.pi
        snaps, _ = allen_cahn(
            N=N, L=L, T_end=0.5, eps=0.15, dt=0.01, save_every=10, seed=42,
        )
        # With random IC small amplitude, solution stays bounded
        assert np.max(np.abs(snaps)) < 5.0, "Solution blew up"


# ======================================================================
# Swift-Hohenberg: f_t = r f − (1+Δ)² f − f³
# ======================================================================


class TestSwiftHohenberg:
    def test_runs_and_produces_output(self):
        N, L = 16, 8 * np.pi
        snaps, times = swift_hohenberg(
            N=N, L=L, T_end=0.01, r=0.3, dt=0.01, save_every=1, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))

    def test_r_negative_returns_to_zero(self):
        """For r < 0, the trivial solution f=0 is stable; small perturbations
        should decay."""
        N, L = 32, 8 * np.pi
        snaps, _ = swift_hohenberg(
            N=N, L=L, T_end=0.5, r=-0.2, dt=0.01, save_every=10, seed=42,
        )
        # L2 norm should decrease from IC (which is small random noise)
        norm0 = np.linalg.norm(snaps[0])
        norm_end = np.linalg.norm(snaps[-1])
        assert norm_end < norm0 * 1.1, (
            f"Solution grew when r<0: {norm0:.4f} → {norm_end:.4f}"
        )

    def test_bounded_solution(self):
        """Swift-Hohenberg solutions should remain bounded."""
        N, L = 32, 8 * np.pi
        snaps, _ = swift_hohenberg(
            N=N, L=L, T_end=0.5, r=0.3, dt=0.01, save_every=10, seed=42,
        )
        assert np.max(np.abs(snaps)) < 10.0, "Solution blew up"


# ======================================================================
# Kuramoto-Sivashinsky: f_t = −Δf − Δ²f − ‖∇f‖²
# ======================================================================


class TestKuramotoSivashinsky:
    def test_runs_and_produces_output(self):
        N, L = 16, 4 * np.pi
        snaps, times = kuramoto_sivashinsky(
            N=N, L=L, T_end=0.05, dt=0.05, save_every=1, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))

    def test_mean_is_not_conserved(self):
        """The spatial mean is NOT conserved for this 2-D (gradient-squared)
        KS form, unlike the classical 1-D convective form f_t = -f f_x - ...
        where -f f_x = -1/2 (f^2)_x is a total derivative and integrates to
        zero exactly. Here, averaging f_t = -1/2|grad f|^2 - Delta f - Delta^2 f
        over the periodic domain kills the Delta/Delta^2 terms (pure
        derivatives integrate to zero) but NOT the -1/2|grad f|^2 term, since
        |grad f|^2 >= 0 always (zero only for constant f). So
            d(mean f)/dt = -1/2 * mean(|grad f|^2) <= 0
        i.e. the mean should decrease monotonically (up to numerical noise),
        never conserved. Verified analytically and numerically (see the
        KS integrator audit in test_ks_integrator.py)."""
        N, L = 32, 4 * np.pi
        snaps, times = kuramoto_sivashinsky(
            N=N, L=L, T_end=5.0, dt=0.02, save_every=5, seed=42,
        )
        means = np.array([np.mean(s) for s in snaps])

        # Should decrease (allow tiny numerical noise, not growth).
        assert means[-1] < means[0] - 1e-6, (
            f"mean should decrease: {means[0]:.6f} -> {means[-1]:.6f}"
        )
        assert np.all(np.diff(means) < 1e-8), "mean should be non-increasing"

        # Quantitatively match the analytic rate d(mean)/dt = -1/2 mean(|grad f|^2)
        # via a central finite difference at an interior snapshot.
        k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
        KX, KY = np.meshgrid(k, k, indexing="ij")
        i = len(snaps) // 2
        fhat = np.fft.fft2(snaps[i])
        fx = np.real(np.fft.ifft2(1j * KX * fhat))
        fy = np.real(np.fft.ifft2(1j * KY * fhat))
        predicted_rate = -0.5 * np.mean(fx**2 + fy**2)
        observed_rate = (means[i + 1] - means[i - 1]) / (times[i + 1] - times[i - 1])
        assert observed_rate == pytest.approx(predicted_rate, rel=0.05)

    def test_bounded_solution(self):
        N, L = 32, 4 * np.pi
        snaps, _ = kuramoto_sivashinsky(
            N=N, L=L, T_end=0.5, dt=0.02, save_every=5, seed=42,
        )
        assert np.max(np.abs(snaps)) < 20.0, "Solution blew up"


# ======================================================================
# Cahn-Hilliard: f_t = −Δ²f + Δ(f³ − f)
# ======================================================================


class TestCahnHilliard:
    def test_runs_and_produces_output(self):
        N, L = 16, 2 * np.pi
        snaps, times = cahn_hilliard(
            N=N, L=L, T_end=0.01, eps=0.2, dt=0.005, save_every=2, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))

    def test_mass_conserved(self):
        """Cahn-Hilliard conserves total mass ∫f dx."""
        N, L = 32, 2 * np.pi
        snaps, _ = cahn_hilliard(
            N=N, L=L, T_end=0.5, eps=0.2, dt=0.01, save_every=10, seed=42,
        )
        masses = [np.sum(s) for s in snaps]
        rel_drift = (max(masses) - min(masses)) / abs(masses[0])
        assert rel_drift < 1e-10, f"Mass drift: {rel_drift:.2e}"

    def test_bounded_solution(self):
        N, L = 32, 2 * np.pi
        snaps, _ = cahn_hilliard(
            N=N, L=L, T_end=0.5, eps=0.2, dt=0.01, save_every=10, seed=42,
        )
        assert np.max(np.abs(snaps)) < 5.0, "Solution blew up"


# ======================================================================
# Porous Medium: f_t = ∇·(f²∇f)
# ======================================================================


class TestPorousMedium:
    def test_runs_and_produces_output(self):
        N, L = 16, 2 * np.pi
        snaps, times = porous_medium(
            N=N, L=L, T_end=0.01, dt=0.005, save_every=2, seed=0,
        )
        assert snaps.ndim == 3
        assert snaps.shape[1] == N
        assert len(times) == snaps.shape[0]
        assert np.all(np.isfinite(snaps))
        # PME should preserve non-negativity
        assert np.all(snaps >= 0)

    def test_mass_conserved(self):
        """PME with periodic BC conserves total mass."""
        N, L = 32, 2 * np.pi
        snaps, _ = porous_medium(
            N=N, L=L, T_end=0.5, dt=0.01, save_every=10, seed=42,
        )
        masses = [np.sum(s) for s in snaps]
        rel_drift = (max(masses) - min(masses)) / abs(masses[0])
        assert rel_drift < 0.01, f"Mass drift: {rel_drift:.2e}"

    def test_non_negativity_preserved(self):
        """PME solution should stay non-negative for all time."""
        N, L = 32, 2 * np.pi
        snaps, _ = porous_medium(
            N=N, L=L, T_end=0.5, dt=0.01, save_every=10, seed=42,
        )
        assert np.all(snaps >= -1e-12), "PME produced negative values"

    def test_maximum_principle(self):
        """The maximum of the PME solution should not increase."""
        N, L = 32, 2 * np.pi
        snaps, _ = porous_medium(
            N=N, L=L, T_end=0.5, dt=0.01, save_every=10, seed=42,
        )
        max0 = np.max(snaps[0])
        for s in snaps[1:]:
            assert np.max(s) <= max0 * 1.01, (
                "PME maximum increased (violates maximum principle)"
            )
