"""Tests for utils.py: print_pde and rotation_90."""

import os
import sys

import numpy as np
import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.utils import print_pde, rotation_90


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
