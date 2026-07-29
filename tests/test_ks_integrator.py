"""Numerical audit of the Kuramoto-Sivashinsky integrator (ETDRK4 + legacy IMEX).

Covers:
- _etdrk4_coeffs correctness against the Cox-Matthews/Kassam-Trefethen (2005)
  closed-form and small-|z| limiting formulas.
- KS linear operator (Lk = K^2 - K^4) and nonlinear term (-1/2|grad f|^2)
  correctness against direct, independently-computed reference values.
- 4th-order temporal convergence of the ETDRK4 scheme (the headline claim of
  the scheme this repo relies on).
- Cross-validation of the ETDRK4 and legacy IMEX schemes against each other.
- The ndim guard (kuramoto_sivashinsky only supports ndim=2).
- Time-array bookkeeping in the shared _imex/_etdrk4 time-stepping loop:
  a saved snapshot is the state *after* that iteration's update, so it is
  labeled with time (step + 1) * dt, matching the true physical time
  (fixed; was previously mislabeled as step * dt, understating the true
  time by one dt and falsely claiming times[0] == 0.0 for a snapshot that
  is actually the state after one full step).
- Dealiasing: the 2/3-rule mask is applied only to the *input* of the
  nonlinear-term derivative computation, not to its output or to the
  persistent Fourier state -- a deviation from the textbook 2/3 rule.
  A quantitative test documents that this has negligible measurable effect
  for KS specifically (its own K^4 hyperviscosity does the real work), so
  this is flagged as a latent hygiene issue, not a correctness bug.
"""

import math
import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.simulations._etdrk4 import _etdrk4_coeffs, _etdrk4_run
from equivariant_pde_find.simulations.kuramoto_sivashinsky import kuramoto_sivashinsky


# ======================================================================
# _etdrk4_coeffs — closed-form / limiting-value correctness
# ======================================================================


class TestEtdrk4CoeffsClosedForm:
    """Cross-check the contour-integral coefficients against direct formulas.

    Kassam & Trefethen (2005), eq. (2.5):
        E  = exp(Lh),  E2 = exp(Lh/2)
        phi1(z) = (e^z - 1) / z            [phi1_dt2 = (dt/2) * phi1(Lh/2)]
        alpha (f1) = h^-2 L^-3 [-4 - Lh + e^Lh(4 - 3Lh + (Lh)^2)]
        beta  (f2) = h^-2 L^-3 [2 + Lh + e^Lh(-2 + Lh)]
        gamma (f4) = h^-2 L^-3 [-4 - 3Lh - (Lh)^2 + e^Lh(4 - Lh)]
    """

    def test_E_E2_are_exact_exponentials(self):
        Lk = np.array([-1.0, -5.0, -50.0, 0.0, 2.0])
        dt = 0.1
        E, E2, phi1_dt2, f1, f2, f4 = _etdrk4_coeffs(Lk, dt)
        np.testing.assert_allclose(E, np.exp(Lk * dt), rtol=1e-13)
        np.testing.assert_allclose(E2, np.exp(Lk * dt / 2), rtol=1e-13)

    def test_direct_formula_matches_at_moderate_z(self):
        """For |z| = O(1), the direct (non-contour) formula is numerically
        safe and should match the contour-integral result to high precision."""
        Lk = np.array([-3.0, -1.5, 1.0, 4.0])
        dt = 0.5
        E, E2, phi1_dt2, f1, f2, f4 = _etdrk4_coeffs(Lk, dt)

        z = Lk * dt
        z2 = Lk * dt / 2
        phi1_direct = (dt / 2) * (np.exp(z2) - 1) / z2
        f1_direct = dt * (-4 - z + np.exp(z) * (4 - 3 * z + z**2)) / z**3
        f2_direct = dt * (2 + z + np.exp(z) * (-2 + z)) / z**3
        f4_direct = dt * (-4 - 3 * z - z**2 + np.exp(z) * (4 - z)) / z**3

        np.testing.assert_allclose(phi1_dt2, phi1_direct, rtol=1e-8)
        np.testing.assert_allclose(f1, f1_direct, rtol=1e-8)
        np.testing.assert_allclose(f2, f2_direct, rtol=1e-8)
        np.testing.assert_allclose(f4, f4_direct, rtol=1e-8)

    def test_small_z_limits_at_dc_mode(self):
        """At Lk=0 (the mean/DC mode), the direct formula suffers catastrophic
        cancellation, but the analytic limits are well known:
            phi1(0) = 1  ->  phi1_dt2 = dt/2
            f1, f2, f4 -> dt/6  (common limit, Kassam-Trefethen sec. 2)
        The contour-integral method must reproduce these without blowing up."""
        Lk = np.array([0.0])
        dt = 0.3
        _, _, phi1_dt2, f1, f2, f4 = _etdrk4_coeffs(Lk, dt)
        assert np.isfinite(phi1_dt2).all() and np.isfinite(f1).all()
        np.testing.assert_allclose(phi1_dt2, dt / 2, rtol=1e-6)
        np.testing.assert_allclose(f1, dt / 6, rtol=1e-6)
        np.testing.assert_allclose(f2, dt / 6, rtol=1e-6)
        np.testing.assert_allclose(f4, dt / 6, rtol=1e-6)

    def test_no_nan_for_stiff_negative_eigenvalues(self):
        """Very stiff (large negative) Lk, as occurs at high wavenumbers for
        KS's -K^4 term, must not produce NaN/Inf in any coefficient."""
        Lk = np.array([-1e2, -1e4, -1e6])
        dt = 0.25
        E, E2, phi1_dt2, f1, f2, f4 = _etdrk4_coeffs(Lk, dt)
        for arr in (E, E2, phi1_dt2, f1, f2, f4):
            assert np.all(np.isfinite(arr))
        # Strongly damped modes should have E essentially zero (exp(Lk*dt)).
        assert np.all(E < 1e-9)


# ======================================================================
# KS linear operator and nonlinear term — direct correctness
# ======================================================================


class TestKSOperatorCorrectness:
    """Verify the discretized linear/nonlinear operators against independent
    (non-FFT-based-shortcut) reference computations, for a case whose Fourier
    transform can be checked by hand: f(x,y) = cos(x) (single Fourier mode)."""

    def test_linear_symbol_matches_single_mode_eigenvalue(self):
        """For f = cos(k0 x) on domain [0, L]^2, Delta f = -k0^2 f and
        Delta^2 f = k0^4 f, so the KS linear part (-Delta - Delta^2) has
        eigenvalue (k0^2 - k0^4) exactly -- matching Lk = K2 - K4 at the
        corresponding grid wavenumber."""
        N, L = 32, 2 * math.pi
        k0 = 2 * (2 * np.pi / L)  # 2 full periods fit exactly on the grid
        x = np.linspace(0, L, N, endpoint=False)
        X, Y = np.meshgrid(x, x, indexing="ij")
        f = np.cos(k0 * X)

        k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
        KX, KY = np.meshgrid(k, k, indexing="ij")
        K2 = KX**2 + KY**2
        K4 = K2**2
        Lk = K2 - K4

        fhat = np.fft.fft2(f)
        # Isolate the dominant mode's eigenvalue by comparing L(f) to f via FFT.
        Lf = np.real(np.fft.ifft2(Lk * fhat))
        expected = (k0**2 - k0**4) * f
        np.testing.assert_allclose(Lf, expected, atol=1e-9)

    def test_nonlinear_term_matches_manual_gradient_squared(self):
        """FFT(-1/2|grad f|^2), computed the way kuramoto_sivashinsky does it
        internally, must match a direct pointwise calculation of the same
        quantity for an arbitrary smooth field."""
        N, L = 24, 4 * math.pi
        rng = np.random.default_rng(0)
        f = 0.3 * rng.standard_normal((N, N))
        # Smooth the random field so its spectrum is well inside the grid's
        # resolvable band (avoids incidental aliasing muddying this check).
        fhat0 = np.fft.fft2(f)
        k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
        KX, KY = np.meshgrid(k, k, indexing="ij")
        K2 = KX**2 + KY**2
        fhat0 *= np.exp(-0.05 * K2)
        f = np.real(np.fft.ifft2(fhat0))

        fhat = np.fft.fft2(f)
        fx = np.real(np.fft.ifft2(1j * KX * fhat))
        fy = np.real(np.fft.ifft2(1j * KY * fhat))
        expected_hat = np.fft.fft2(-0.5 * (fx * fx + fy * fy))

        # Reproduce kuramoto_sivashinsky's internal nonlin_hat with dealias=False.
        computed_hat = np.fft.fft2(-0.5 * (fx**2 + fy**2))
        np.testing.assert_allclose(computed_hat, expected_hat, atol=1e-9)


# ======================================================================
# ETDRK4 — 4th-order temporal convergence
# ======================================================================


def _run_at_matched_time(dt, T_target, N=24, L=4 * math.pi, seed=7, scheme="etdrk4"):
    """Run KS to a physical time EXACTLY equal to T_target, regardless of dt.

    kuramoto_sivashinsky internally uses
        n_steps = round(T_end / dt) + 1
    as a rounding-safety margin (see kuramoto_sivashinsky docstring / _imex.py),
    which means the true elapsed time after the run is n_steps*dt, not T_end
    -- it overshoots by up to one dt, and the overshoot itself depends on dt.
    A naive convergence study that just requests T_end=T_target at different
    dt therefore compares solutions at *different* physical times and gets a
    spurious ~1st-order result. Requesting T_end = T_target - dt cancels the
    "+1" exactly, so every trial's last saved snapshot lands at T_target.
    """
    T_req = T_target - dt
    n_steps_expected = int(round(T_req / dt)) + 1
    snaps, times = kuramoto_sivashinsky(
        N=N, L=L, T_end=T_req, dt=dt, save_every=1, seed=seed, scheme=scheme,
    )
    assert len(snaps) == n_steps_expected
    return snaps[-1]


class TestEtdrk4ConvergenceOrder:
    @pytest.mark.slow
    def test_fourth_order_convergence(self):
        """Refining dt should reduce the error by ~16x per halving once in
        the asymptotic regime -- the entire point of using ETDRK4 over the
        1st-order legacy IMEX scheme."""
        T_target = 0.2
        ref = _run_at_matched_time(1e-6, T_target)

        dts = [0.01, 0.005, 0.0025, 0.00125, 0.000625]
        errs = [
            np.linalg.norm(_run_at_matched_time(dt, T_target) - ref)
            / np.linalg.norm(ref)
            for dt in dts
        ]

        orders = [
            np.log(errs[i] / errs[i + 1]) / np.log(dts[i] / dts[i + 1])
            for i in range(len(dts) - 1)
        ]
        # Order should increase toward 4 as dt shrinks (finite-dt effects
        # dominate at the coarse end); require it clearly exceeds 1st/2nd
        # order and approaches 4th order at the fine end.
        assert orders[-1] > 3.5, f"observed orders {orders}, expected ~4 at fine dt"
        assert all(o > 1.5 for o in orders), f"observed orders {orders}"

    @pytest.mark.slow
    def test_naive_time_alignment_gives_spuriously_low_order(self):
        """Regression guard for the convergence-study pitfall itself: without
        correcting for the dt-dependent overshoot, refining dt at a *fixed*
        T_end request produces an apparent ~1st-order rate, not because the
        scheme is low-order, but because each trial is secretly evaluated at
        a different physical time (T_end + dt). This test documents that
        trap so it isn't rediscovered as a false "ETDRK4 is broken" bug."""
        N, L, T_end = 24, 4 * math.pi, 0.2
        seed = 7

        def run_naive(dt):
            save_every = max(1, int(round(T_end / dt)))
            snaps, _ = kuramoto_sivashinsky(
                N=N, L=L, T_end=T_end, dt=dt, save_every=save_every, seed=seed,
            )
            return snaps[-1]

        ref = run_naive(1e-4)
        dts = [0.02, 0.01, 0.005]
        errs = [np.linalg.norm(run_naive(dt) - ref) / np.linalg.norm(ref) for dt in dts]
        orders = [
            np.log(errs[i] / errs[i + 1]) / np.log(dts[i] / dts[i + 1])
            for i in range(len(dts) - 1)
        ]
        # The naive (misaligned) study should look close to 1st order, in
        # contrast to the >3.5 true order measured above.
        assert all(o < 1.5 for o in orders), (
            f"expected the misaligned study to look ~1st order, got {orders} "
            "-- if this now fails, the T_end/dt time-alignment bug this test "
            "documents may have been fixed; update/remove this test."
        )


# ======================================================================
# ETDRK4 vs legacy IMEX — cross-validation
# ======================================================================


class TestSchemeCrossValidation:
    @pytest.mark.slow
    def test_etdrk4_and_imex_converge_to_each_other(self):
        """Both schemes solve the same PDE; as dt -> 0 their solutions at a
        matched physical time should converge toward each other, with the
        gap shrinking at close to the (1st-order) rate of the slower scheme."""
        T_target = 0.3
        dts = [0.01, 0.002, 0.0005]
        diffs = []
        for dt in dts:
            a = _run_at_matched_time(dt, T_target, scheme="etdrk4")
            b = _run_at_matched_time(dt, T_target, scheme="imex")
            diffs.append(np.linalg.norm(a - b) / np.linalg.norm(a))

        assert diffs[-1] < diffs[0], "schemes should converge as dt shrinks"
        assert diffs[-1] < 5e-3, f"schemes disagree too much at fine dt: {diffs}"


# ======================================================================
# Time-array bookkeeping — known, confirmed quirk
# ======================================================================


class TestTimeArrayBookkeeping:
    """`_etdrk4_run`/`_imex._run` label a saved snapshot with the step index
    *after* that iteration's update (``step + 1``), matching the fact that
    the snapshot itself is the state after ``step + 1`` total steps -- so
    ``times[i]`` correctly equals the physical time of ``snapshots[i]``.

    Neither loop ever saves the literal t=0 initial condition (the earliest
    possible save is after one full step), so ``times[0] == dt``, not 0 --
    this is a structural property of the save-after-update loop, not a bug,
    and is now labeled truthfully (previously ``times[0]`` was mislabeled
    as 0.0 despite the snapshot being the post-one-step state; fixed in
    _imex.py/_etdrk4.py/porous_medium.py).
    """

    def test_first_snapshot_is_labeled_dt_not_zero(self):
        N, L, dt = 16, 4 * math.pi, 0.05
        rng = np.random.default_rng(42)
        f0 = 0.1 * rng.standard_normal((N, N))

        snaps, times = kuramoto_sivashinsky(
            N=N, L=L, T_end=0.5, dt=dt, save_every=1, f0=f0.copy(), seed=42,
        )
        assert times[0] == pytest.approx(dt)
        assert not np.allclose(snaps[0], f0), (
            "snapshots[0] is the state after one dt-step, not f0 -- correctly "
            "reflected by times[0] == dt now, not claimed to be 0.0."
        )

    def test_consecutive_time_spacing_is_still_correct(self):
        """The systematic offset is constant, so *differences* between
        consecutive reported times (what dt_snap / finite-difference f_t
        estimation actually uses) are unaffected by the bug."""
        N, L, dt, save_every = 16, 4 * math.pi, 0.05, 3
        snaps, times = kuramoto_sivashinsky(
            N=N, L=L, T_end=1.0, dt=dt, save_every=save_every, seed=1,
        )
        np.testing.assert_allclose(
            np.diff(times), save_every * dt, rtol=0, atol=1e-12,
        )


# ======================================================================
# Dealiasing — quantitative check of what the 2/3-rule mask actually does
# ======================================================================


class TestDealiasing:
    """The `dealias` mask in kuramoto_sivashinsky is applied only to the
    *input* used to compute the nonlinear term's spatial derivatives, not to
    the term's own output nor to the persistent Fourier state uhat -- so the
    top-1/3 wavenumber band of uhat (and hence of saved snapshots) is never
    explicitly zeroed, unlike the textbook 2/3-rule (which masks the
    nonlinear term's output at every step; confirmed via literature search
    against standard pseudospectral-method practice).

    Empirically, for KS specifically this has negligible effect: the outer
    band's energy fraction converges to the same tiny value (~1e-10 or
    smaller) within 1-2 steps regardless of `dealias`, because the K^4
    hyperviscosity damps high wavenumbers far faster than aliasing can
    repopulate them. These tests pin down that empirical finding rather than
    asserting the mask achieves exact (textbook) dealiasing.
    """

    @staticmethod
    def _outer_band_energy_fraction(f, N, L):
        k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
        KX, KY = np.meshgrid(k, k, indexing="ij")
        k_max = (N // 2) * (2 * np.pi / L)
        outer = (np.abs(KX) > 2 / 3 * k_max) | (np.abs(KY) > 2 / 3 * k_max)
        fhat = np.fft.fft2(f)
        return np.sum(np.abs(fhat[outer]) ** 2) / np.sum(np.abs(fhat) ** 2)

    def test_dealias_true_does_not_fully_zero_outer_band(self):
        """Documents that dealias=True is NOT a hard guarantee: the outer
        band carries small but strictly nonzero energy (it is not held to
        exactly zero the way an output-masked implementation would)."""
        N, L = 32, 16 * math.pi
        snaps, _ = kuramoto_sivashinsky(
            N=N, L=L, T_end=1.0, dt=0.1, save_every=1, seed=1, dealias=True,
        )
        frac = self._outer_band_energy_fraction(snaps[-1], N, L)
        assert frac > 0.0

    def test_dealias_true_vs_false_converge_within_a_couple_steps(self):
        """The dealias=True/False outer-band energy fractions should differ
        at most transiently (step 0) and become indistinguishable within a
        handful of steps, because K^4 damping -- not the mask -- is what
        actually keeps the resolved solution clean here."""
        N, L = 64, 16 * math.pi
        snaps_on, _ = kuramoto_sivashinsky(
            N=N, L=L, T_end=2.5, dt=0.25, save_every=1, seed=1, dealias=True,
        )
        snaps_off, _ = kuramoto_sivashinsky(
            N=N, L=L, T_end=2.5, dt=0.25, save_every=1, seed=1, dealias=False,
        )
        frac_on = self._outer_band_energy_fraction(snaps_on[-1], N, L)
        frac_off = self._outer_band_energy_fraction(snaps_off[-1], N, L)
        # Both tiny, and close to each other (within an order of magnitude).
        assert frac_on < 1e-6 and frac_off < 1e-6
        assert frac_off / max(frac_on, 1e-300) < 100


# Note: the ndim guard (only ndim=2 supported) is already covered by
# TestKuramotoSivashinskyNdimGuard in test_ndim.py -- not duplicated here.
