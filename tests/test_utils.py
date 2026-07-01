"""Tests for utils.py: print_pde, rotation_90, and anisotropic noise."""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.utils import (
    print_pde, rotation_90,
    make_spatial_grid, build_noise_profile, normalize_profile,
    apply_anisotropic_noise,
)


class TestPrintPde:
    def test_empty_coefficients(self, capsys):
        w = np.zeros(5)
        descs = ["f", "Δf", "‖∇f‖²", "f³", "Δ²f"]
        print_pde(w, descs)
        captured = capsys.readouterr()
        assert "f_t = 0" in captured.out

    def test_single_term(self, capsys):
        w = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
        descs = ["f", "Δf", "‖∇f‖²", "f³", "Δ²f"]
        print_pde(w, descs)
        captured = capsys.readouterr()
        assert "f_t" in captured.out
        assert "+1.00000" in captured.out
        assert "f" in captured.out

    def test_multiple_terms(self, capsys):
        w = np.array([1.0, 0.0, -1.0, -1.0, 0.0])
        descs = ["f", "Δf", "‖∇f‖²", "f³", "Δ²f"]
        print_pde(w, descs)
        captured = capsys.readouterr()
        assert "+" in captured.out
        assert "-1.00000" in captured.out

    def test_custom_lhs(self, capsys):
        w = np.array([1.0, 0.0, 0.0, 0.0, 0.0])
        descs = ["f", "Δf", "‖∇f‖²", "f³", "Δ²f"]
        print_pde(w, descs, lhs="u_t")
        captured = capsys.readouterr()
        assert "u_t" in captured.out


class TestRotation90:
    def test_square_field(self):
        N = 32
        f = np.random.default_rng(0).standard_normal((N, N))
        f_rot = rotation_90(f)

        # Shape preserved
        assert f_rot.shape == (N, N)

        # Rotating 4 times should give back the original
        f_rot4 = rotation_90(rotation_90(rotation_90(f_rot)))
        np.testing.assert_allclose(f_rot4, f, atol=1e-14)

    def test_rotation_preserves_norm(self):
        N = 32
        f = np.random.default_rng(1).standard_normal((N, N))
        f_rot = rotation_90(f)
        assert np.linalg.norm(f_rot) == pytest.approx(np.linalg.norm(f))

    def test_constant_field(self):
        N = 32
        f = np.ones((N, N))
        f_rot = rotation_90(f)
        np.testing.assert_allclose(f_rot, f)


# ── anisotropic noise ───────────────────────────────────────────────────────


class TestSpatialGrid:
    def test_shape(self):
        X, Y = make_spatial_grid(64, 8 * np.pi, 8 * np.pi)
        assert X.shape == (64, 64)
        assert Y.shape == (64, 64)

    def test_centering(self):
        X, Y = make_spatial_grid(64, 10.0, 10.0)
        assert X.min() == pytest.approx(-5.0)
        assert X.max() < 5.0  # endpoint=False
        assert Y.min() == pytest.approx(-5.0)

    def test_ij_indexing(self):
        """With indexing='ij', X varies along axis 0 (rows), Y along axis 1
        (columns).  For fixed i, X[i, :] is constant; for fixed j, Y[:, j] is
        constant."""
        X, Y = make_spatial_grid(32, 4.0, 4.0)
        # X[i, j] should equal X[i, 0] for all j
        for i in range(32):
            np.testing.assert_allclose(X[i, :], X[i, 0])
        # Y[i, j] should equal Y[0, j] for all i
        for j in range(32):
            np.testing.assert_allclose(Y[:, j], Y[0, j])


class TestNoiseProfile:
    def test_radial_symmetry(self):
        """Radial profile depends only on distance from origin.

        On a periodic grid with endpoint=False, the point symmetric to (i, j)
        across x=0 is ((N-i) % N, j), and across y=0 is (i, (N-j) % N).
        The diagonal swap uses the same reflection logic."""
        N = 32
        X, Y = make_spatial_grid(N, 8.0, 8.0)
        P = build_noise_profile(X, Y, direction="radial")
        # Reflection across x = 0
        for i in [1, 5, 16]:
            i_sym = (N - i) % N
            np.testing.assert_allclose(P[i, :], P[i_sym, :], rtol=1e-14)
        # Reflection across y = 0
        for j in [1, 5, 16]:
            j_sym = (N - j) % N
            np.testing.assert_allclose(P[:, j], P[:, j_sym], rtol=1e-14)
        # Diagonal swap: r(x, y) = r(y, x), but only for points with the
        # same distance — test a few representative pairs
        np.testing.assert_allclose(P[8, 16], P[16, 8], rtol=1e-14)

    def test_axial_depends_only_on_abs_x(self):
        """Axial profile P(x, y) = P(-x, y) and P(x, y₁) = P(x, y₂).

        On a periodic grid with endpoint=False, the x-reflection maps
        index i to (N-i) % N."""
        N = 32
        X, Y = make_spatial_grid(N, 8.0, 8.0)
        P = build_noise_profile(X, Y, direction="axial")
        # x-reflection
        np.testing.assert_allclose(P[1, :], P[(N - 1) % N, :], rtol=1e-14)
        np.testing.assert_allclose(P[5, :], P[(N - 5) % N, :], rtol=1e-14)
        # y-invariance: P[i, j₁] = P[i, j₂] for all j₁, j₂
        for i in [0, 8, 16, 24]:
            np.testing.assert_allclose(P[i, :], P[i, 0], rtol=1e-14)
        # But P(x, y) ≠ P(y, x) in general — spot-check
        assert P[8, 16] != pytest.approx(P[16, 8])

    def test_gaussian_decay_range(self):
        X, Y = make_spatial_grid(32, 8.0, 8.0)
        P = build_noise_profile(X, Y, profile_type="gaussian_decay")
        assert P.max() == pytest.approx(1.0)   # at origin, r=0
        assert P.min() > 0.0                    # never exactly zero
        assert P.min() < P.max()                # decays

    def test_gaussian_growth_range(self):
        X, Y = make_spatial_grid(32, 8.0, 8.0)
        P = build_noise_profile(
            X, Y, profile_type="gaussian_growth", amplitude=2.0,
        )
        assert P.min() == pytest.approx(1.0)    # at origin
        assert P.max() <= 3.0                    # 1 + A at large r
        assert P.max() > P.min()                 # grows outward

    def test_defaults_produce_valid_profile(self):
        X, Y = make_spatial_grid(64, 4 * np.pi, 4 * np.pi)
        for ptype in ["gaussian_decay", "gaussian_growth"]:
            for direction in ["radial", "axial"]:
                P = build_noise_profile(
                    X, Y, profile_type=ptype, direction=direction,
                )
                assert P.shape == (64, 64)
                assert np.all(P >= 0.0)
                assert not np.any(np.isnan(P))

    def test_unknown_direction_raises(self):
        X, Y = make_spatial_grid(8, 1.0, 1.0)
        with pytest.raises(ValueError):
            build_noise_profile(X, Y, direction="diagonal")

    def test_unknown_profile_type_raises(self):
        X, Y = make_spatial_grid(8, 1.0, 1.0)
        with pytest.raises(ValueError):
            build_noise_profile(X, Y, profile_type="exponential")


class TestNormalizeProfile:
    def test_normalization(self):
        X, Y = make_spatial_grid(64, 8.0, 8.0)
        P = build_noise_profile(X, Y, profile_type="gaussian_growth",
                                direction="radial", amplitude=2.0)
        Pn = normalize_profile(P)
        assert np.mean(Pn ** 2) == pytest.approx(1.0, rel=1e-12)

    def test_normalization_independent_of_scale(self):
        X, Y = make_spatial_grid(32, 1.0, 1.0)
        P1 = build_noise_profile(X, Y, profile_type="gaussian_decay",
                                 direction="radial", length_scale=0.25)
        P2 = build_noise_profile(X, Y, profile_type="gaussian_decay",
                                 direction="radial", length_scale=1.0)
        # Different raw profiles should both normalise to unit RMS
        assert np.mean(normalize_profile(P1) ** 2) == pytest.approx(1.0)
        assert np.mean(normalize_profile(P2) ** 2) == pytest.approx(1.0)

    def test_normalized_profile_nonnegative(self):
        X, Y = make_spatial_grid(32, 4.0, 4.0)
        for ptype in ["gaussian_decay", "gaussian_growth"]:
            P = build_noise_profile(X, Y, profile_type=ptype)
            Pn = normalize_profile(P)
            assert np.all(Pn >= 0.0)


class TestApplyAnisotropicNoise:
    def test_shape_preserved(self):
        snapshots = np.random.default_rng(0).standard_normal((10, 32, 32))
        X, Y = make_spatial_grid(32, 4.0, 4.0)
        P = normalize_profile(
            build_noise_profile(X, Y, direction="radial"),
        )
        rng = np.random.default_rng(42)
        noisy = apply_anisotropic_noise(snapshots, 0.02, P, rng)
        assert noisy.shape == snapshots.shape

    def test_noise_power_matches_isotropic(self):
        """Anisotropic noise with normalised profile should have same expected
        Frobenius norm as isotropic noise at the same sigma."""
        snapshots = np.random.default_rng(1).standard_normal((100, 16, 16))
        X, Y = make_spatial_grid(16, 2.0, 2.0)
        P = normalize_profile(
            build_noise_profile(X, Y, direction="axial",
                                profile_type="gaussian_growth"),
        )
        obs_noise = 0.05
        sigma = obs_noise * snapshots.std()

        # Isotropic noise power (theoretical)
        expected_var = sigma ** 2  # per element
        expected_norm_sq = expected_var * snapshots.size

        # Anisotropic noise — average over many draws
        rng = np.random.default_rng(42)
        norms_sq = []
        for _ in range(200):
            noisy = apply_anisotropic_noise(snapshots, obs_noise, P, rng)
            noise_field = noisy - snapshots
            norms_sq.append(np.sum(noise_field ** 2))
        mean_norm_sq = np.mean(norms_sq)

        # Should match within ~1% (200 draws)
        assert mean_norm_sq == pytest.approx(expected_norm_sq, rel=0.02)

    def test_zero_noise_returns_clean(self):
        snapshots = np.random.default_rng(2).standard_normal((5, 16, 16))
        X, Y = make_spatial_grid(16, 2.0, 2.0)
        P = normalize_profile(build_noise_profile(X, Y))
        rng = np.random.default_rng(42)
        noisy = apply_anisotropic_noise(snapshots, 0.0, P, rng)
        np.testing.assert_allclose(noisy, snapshots, atol=1e-14)

    def test_reproducibility(self):
        """Same seed → same noise."""
        snapshots = np.random.default_rng(3).standard_normal((5, 16, 16))
        X, Y = make_spatial_grid(16, 2.0, 2.0)
        P = normalize_profile(build_noise_profile(X, Y))
        rng1 = np.random.default_rng(42)
        rng2 = np.random.default_rng(42)
        noisy1 = apply_anisotropic_noise(snapshots, 0.02, P, rng1)
        noisy2 = apply_anisotropic_noise(snapshots, 0.02, P, rng2)
        np.testing.assert_allclose(noisy1, noisy2, atol=1e-14)
