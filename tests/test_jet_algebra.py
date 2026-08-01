"""Tests for ``jet_algebra.py`` — parse_jet_monomial (reverse parser)."""

import os
import sys

import pytest

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, CODE_DIR)

from equivariant_pde_find.jet_algebra import parse_jet_monomial


# ---------------------------------------------------------------------------
# Basic parsing
# ---------------------------------------------------------------------------


class TestParseJetMonomial:
    """Round-trip tests: parse_jet_monomial should be the inverse of
    _describe_exponents / describe_monomial."""

    def test_constant(self):
        assert parse_jet_monomial("1", 2) == {}

    def test_f(self):
        assert parse_jet_monomial("f", 2) == {(0, 0): 1}

    def test_f_x(self):
        assert parse_jet_monomial("f_x", 2) == {(1, 0): 1}

    def test_f_y(self):
        assert parse_jet_monomial("f_y", 2) == {(0, 1): 1}

    def test_f_xx(self):
        assert parse_jet_monomial("f_xx", 2) == {(2, 0): 1}

    def test_f_xy(self):
        assert parse_jet_monomial("f_xy", 2) == {(1, 1): 1}

    def test_f_yy(self):
        assert parse_jet_monomial("f_yy", 2) == {(0, 2): 1}

    def test_f_x_squared(self):
        assert parse_jet_monomial("f_x²", 2) == {(1, 0): 2}

    def test_f_y_squared(self):
        assert parse_jet_monomial("f_y²", 2) == {(0, 1): 2}

    def test_f_cubed(self):
        assert parse_jet_monomial("f³", 2) == {(0, 0): 3}

    def test_f_times_f_x(self):
        assert parse_jet_monomial("f·f_x", 2) == {(0, 0): 1, (1, 0): 1}

    def test_f_cubed_times_f_x(self):
        assert parse_jet_monomial("f³·f_x", 2) == {(0, 0): 3, (1, 0): 1}

    def test_f_times_f_xx(self):
        assert parse_jet_monomial("f·f_xx", 2) == {(0, 0): 1, (2, 0): 1}

    def test_f_x_squared_times_f_xx(self):
        assert parse_jet_monomial("f_x²·f_xx", 2) == {(1, 0): 2, (2, 0): 1}

    def test_f_squared_times_f_xx(self):
        assert parse_jet_monomial("f²·f_xx", 2) == {(0, 0): 2, (2, 0): 1}

    def test_f_xx_times_f_yy(self):
        assert parse_jet_monomial("f_xx·f_yy", 2) == {(2, 0): 1, (0, 2): 1}

    def test_f_times_f_x_squared(self):
        assert parse_jet_monomial("f·f_x²", 2) == {(0, 0): 1, (1, 0): 2}

    def test_f_times_f_y_squared(self):
        assert parse_jet_monomial("f·f_y²", 2) == {(0, 0): 1, (0, 1): 2}

    def test_fourth_order(self):
        assert parse_jet_monomial("f_xxxx", 2) == {(4, 0): 1}

    def test_mixed_fourth(self):
        assert parse_jet_monomial("f_xxyy", 2) == {(2, 2): 1}

    # ------------------------------------------------------------------
    # ndim=3
    # ------------------------------------------------------------------

    def test_ndim3_f_z(self):
        assert parse_jet_monomial("f_z", 3) == {(0, 0, 1): 1}

    def test_ndim3_f_z_squared(self):
        assert parse_jet_monomial("f_z²", 3) == {(0, 0, 1): 2}

    def test_ndim3_f_xz(self):
        assert parse_jet_monomial("f_xz", 3) == {(1, 0, 1): 1}

    # ------------------------------------------------------------------
    # Error cases
    # ------------------------------------------------------------------

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            parse_jet_monomial("", 2)

    def test_bad_factor_raises(self):
        with pytest.raises(ValueError):
            parse_jet_monomial("g_x", 2)

    def test_bad_axis_raises(self):
        with pytest.raises(ValueError):
            parse_jet_monomial("f_q", 2)

    def test_axis_out_of_range(self):
        with pytest.raises(ValueError):
            parse_jet_monomial("f_z", 2)


# ---------------------------------------------------------------------------
# Round-trip with the polyjet formatter
# ---------------------------------------------------------------------------


class TestRoundTripWithPolyjet:
    """Verify parse_jet_monomial inverts _describe_exponents."""

    @staticmethod
    def _describe_exponents(exponents, jet_names):
        """Vendored from polyjet.py to avoid import coupling."""
        _SUP = {2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹"}

        def _sup(n):
            return _SUP.get(n, f"^{n}")

        parts = []
        for mi, exp in sorted(exponents.items(), key=lambda kv: (sum(kv[0]), kv[0])):
            if exp == 0:
                continue
            name = jet_names[mi]
            parts.append(name if exp == 1 else f"{name}{_sup(exp)}")
        return "1" if not parts else "·".join(parts)

    def _jet_names(self, ndim):
        """Build {multi_index: name} dict."""
        labels = ["x", "y", "z", "w", "v", "u", "t", "s"]
        names = {}
        # f
        names[(0,) * ndim] = "f"
        # f_x, f_xx, f_xxx, f_xxxx, f_xy, etc.
        for total in range(1, 5):
            for mi in self._indices_of_order(ndim, total):
                parts = "".join(labels[i] * a for i, a in enumerate(mi) if a)
                names[mi] = "f_" + parts
        return names

    @staticmethod
    def _indices_of_order(ndim, total):
        if ndim == 1:
            return [(total,)]
        out = []
        for a in range(total + 1):
            for rest in TestRoundTripWithPolyjet._indices_of_order(ndim - 1, total - a):
                out.append((a,) + rest)
        return out

    def test_round_trip_basic_monomials(self):
        """parse → format → same string for common monomials."""
        ndim = 2
        jet_names = self._jet_names(ndim)
        descs = [
            "1", "f", "f_x", "f_y", "f_xx", "f_xy", "f_yy",
            "f_x²", "f_y²", "f³", "f·f_x", "f³·f_x",
            "f·f_xx", "f_x²·f_xx", "f²·f_xx", "f_yy·f_xx",
            "f·f_x²", "f·f_y²", "f_xxxx", "f_xxyy", "f_yyyy",
        ]
        for desc in descs:
            exponents = parse_jet_monomial(desc, ndim)
            formatted = self._describe_exponents(exponents, jet_names)
            assert formatted == desc, f"Round-trip failed: {desc!r} → {exponents} → {formatted!r}"

    def test_round_trip_ndim3(self):
        """parse → format → same string for 3-D monomials."""
        ndim = 3
        jet_names = self._jet_names(ndim)
        descs = ["1", "f", "f_z", "f_z²", "f_xz", "f·f_z", "f_xxz"]
        for desc in descs:
            exponents = parse_jet_monomial(desc, ndim)
            formatted = self._describe_exponents(exponents, jet_names)
            assert formatted == desc, f"Round-trip failed: {desc!r} → {exponents} → {formatted!r}"
