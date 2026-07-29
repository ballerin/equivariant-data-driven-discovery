"""Tests for ``equivariant_pde_find.equations``."""

import numpy as np
import pytest

from equivariant_pde_find.equations import (
    CH_EPSILON,
    EQUATION_META,
    _EQ_NAME_ALIASES,
    _pdefind_axis_labels,
    find_op,
    mask_pme,
    normalize_eq_name,
    true_allen_cahn,
    true_allen_cahn_mod_multigraph,
    true_allen_cahn_pdefind,
    true_cahn_hilliard,
    true_cahn_hilliard_mod_multigraph,
    true_cahn_hilliard_pdefind,
    true_ks,
    true_ks_mod_multigraph,
    true_pme,
    true_pme_mod_multigraph,
    true_pme_wsindy,
    true_swift_hohenberg,
    true_swift_hohenberg_mod_multigraph,
)

# ---------------------------------------------------------------------------
# Mock description lists
# ---------------------------------------------------------------------------

# A plausible mock for the original multigraph basis (20 ops).
_MOCK_DESCRIPTIONS_EQUIVARIANT = [
    '1', 'f', 'f²', 'f³',
    'Δf', 'Δ²f', '‖∇f‖²', '‖∇²f‖²',
    'f·Δf', 'f²·Δf', 'f·Δ²f', 'f²·Δ²f',
    '⟨∇Δf,∇f⟩', 'f·⟨∇Δf,∇f⟩ + ‖∇f‖²·Δf',
    'f·‖∇f‖²', 'Δf³',
    '∇f·∇(Δf)', '∇·(f∇f)', '∇·(f²∇f)', '∇·(‖∇f‖²∇f)',
]

# Mock for modified multigraph basis — same length, different col 11/13 → Δf³ / Δ(f·‖∇f‖²)
_MOCK_DESCRIPTIONS_MOD = list(_MOCK_DESCRIPTIONS_EQUIVARIANT)
_MOCK_DESCRIPTIONS_MOD[11] = 'Δf³'
_MOCK_DESCRIPTIONS_MOD[13] = 'Δ(f·‖∇f‖²)'  # placeholder

# Mock for WSINDy basis
_MOCK_DESCRIPTIONS_WSINDY = [
    'f', 'f\u00b2', 'f\u00b3',
    'f_x', 'f_y', 'f_xx', 'f_xy', 'f_yy',
    '(f\u00b2)_x', '(f\u00b2)_y', '(f\u00b2)_xx', '(f\u00b2)_xy', '(f\u00b2)_yy',
    '(f\u00b3)_x', '(f\u00b3)_y', '(f\u00b3)_{xx}', '(f\u00b3)_{xy}', '(f\u00b3)_{yy}',
]


# ---------------------------------------------------------------------------
# Tests — constants
# ---------------------------------------------------------------------------

def test_ch_epsilon():
    assert CH_EPSILON == 0.2


# ---------------------------------------------------------------------------
# Tests — find_op
# ---------------------------------------------------------------------------

def test_find_op_returns_correct_index():
    assert find_op(_MOCK_DESCRIPTIONS_EQUIVARIANT, 'f') == 1
    assert find_op(_MOCK_DESCRIPTIONS_EQUIVARIANT, 'Δ²f') == 5


def test_find_op_raises_keyerror_for_missing():
    with pytest.raises(KeyError, match="not in library"):
        find_op(_MOCK_DESCRIPTIONS_EQUIVARIANT, 'nonexistent')


# ---------------------------------------------------------------------------
# Tests — mask_pme
# ---------------------------------------------------------------------------

def test_mask_pme_filters_near_zero():
    f = np.array([[0.0, 0.03], [0.05, 0.10]])
    mask = mask_pme(f)
    assert mask.shape == (4,)
    assert mask.dtype == bool
    # 0.0 < 0.05, 0.03 < 0.05, 0.05 is NOT > 0.05, 0.10 > 0.05
    expected = [False, False, False, True]
    np.testing.assert_array_equal(mask, expected)


# ---------------------------------------------------------------------------
# Tests — _pdefind_axis_labels
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ndim,expected", [
    (1, ['x']),
    (2, ['x', 'y']),
    (3, ['x', 'y', 'z']),
    (8, ['x', 'y', 'z', 'w', 'v', 'u', 't', 's']),
])
def test_pdefind_axis_labels_valid(ndim, expected):
    assert _pdefind_axis_labels(ndim) == expected


def test_pdefind_axis_labels_raises_for_invalid_ndim():
    with pytest.raises(ValueError, match="ndim must be 1-8"):
        _pdefind_axis_labels(0)
    with pytest.raises(ValueError, match="ndim must be 1-8"):
        _pdefind_axis_labels(9)


# ---------------------------------------------------------------------------
# Tests — truth functions: original multigraph basis
# ---------------------------------------------------------------------------

def test_true_allen_cahn():
    desc = _MOCK_DESCRIPTIONS_EQUIVARIANT
    w = true_allen_cahn(desc, eps=0.15)
    assert len(w) == len(desc)
    assert w[find_op(desc, 'f')] == 1.0
    assert w[find_op(desc, 'f³')] == -1.0
    assert w[find_op(desc, 'Δf')] == pytest.approx(0.0225)
    # All other entries are zero
    nonzero = set(np.nonzero(w)[0])
    assert nonzero == {find_op(desc, 'f'), find_op(desc, 'f³'), find_op(desc, 'Δf')}


def test_true_swift_hohenberg():
    desc = _MOCK_DESCRIPTIONS_EQUIVARIANT
    w = true_swift_hohenberg(desc, r=0.3)
    assert w[find_op(desc, 'f')] == -0.7
    assert w[find_op(desc, 'Δf')] == -2.0
    assert w[find_op(desc, 'Δ²f')] == -1.0
    assert w[find_op(desc, 'f³')] == -1.0
    assert np.count_nonzero(w) == 4


def test_true_ks():
    desc = _MOCK_DESCRIPTIONS_EQUIVARIANT
    w = true_ks(desc)
    assert w[find_op(desc, '‖∇f‖²')] == -0.5
    assert w[find_op(desc, 'Δf')] == -1.0
    assert w[find_op(desc, 'Δ²f')] == -1.0
    assert np.count_nonzero(w) == 3


def test_true_cahn_hilliard():
    desc = _MOCK_DESCRIPTIONS_EQUIVARIANT
    w = true_cahn_hilliard(desc, eps=0.2)
    assert w[find_op(desc, 'Δf')] == -1.0
    assert w[find_op(desc, 'f²·Δf')] == 3.0
    assert w[find_op(desc, 'f·‖∇f‖²')] == 6.0
    assert w[find_op(desc, 'Δ²f')] == pytest.approx(-0.04)
    assert np.count_nonzero(w) == 4


def test_true_pme():
    desc = _MOCK_DESCRIPTIONS_EQUIVARIANT
    w = true_pme(desc)
    assert w[find_op(desc, 'f²·Δf')] == 1.0
    assert w[find_op(desc, 'f·‖∇f‖²')] == 2.0
    assert np.count_nonzero(w) == 2


# ---------------------------------------------------------------------------
# Tests — truth functions: modified multigraph basis
# ---------------------------------------------------------------------------

def test_mod_multigraph_identity_for_ac_sh_ks():
    """For AC, SH, KS the modified truth equals the original truth."""
    desc_orig = _MOCK_DESCRIPTIONS_EQUIVARIANT
    desc_mod = _MOCK_DESCRIPTIONS_MOD
    np.testing.assert_array_equal(
        true_allen_cahn_mod_multigraph(desc_mod),
        true_allen_cahn(desc_orig),
    )
    np.testing.assert_array_equal(
        true_swift_hohenberg_mod_multigraph(desc_mod),
        true_swift_hohenberg(desc_orig),
    )
    np.testing.assert_array_equal(
        true_ks_mod_multigraph(desc_mod),
        true_ks(desc_orig),
    )


def test_true_cahn_hilliard_mod_multigraph():
    desc = _MOCK_DESCRIPTIONS_MOD
    w = true_cahn_hilliard_mod_multigraph(desc, eps=0.2)
    assert w[find_op(desc, 'Δf')] == -1.0
    assert w[find_op(desc, 'Δf³')] == 1.0
    assert w[find_op(desc, 'Δ²f')] == pytest.approx(-0.04)
    # f·‖∇f‖² should be zero in the modified basis
    assert w[find_op(desc, 'f·‖∇f‖²')] == 0.0
    assert np.count_nonzero(w) == 3


def test_true_pme_mod_multigraph():
    desc = _MOCK_DESCRIPTIONS_MOD
    w = true_pme_mod_multigraph(desc)
    assert w[find_op(desc, 'Δf³')] == pytest.approx(1.0 / 3.0)
    assert w[find_op(desc, 'f·‖∇f‖²')] == 0.0
    assert np.count_nonzero(w) == 1


# ---------------------------------------------------------------------------
# Tests — truth functions: PDE-FIND basis
# ---------------------------------------------------------------------------

def test_true_allen_cahn_pdefind_2d():
    result = true_allen_cahn_pdefind(2, eps=0.15)
    assert result["f"] == 1.0
    assert result["f³"] == -1.0
    assert result["f_xx"] == pytest.approx(0.0225)
    assert result["f_yy"] == pytest.approx(0.0225)
    assert len(result) == 4  # f, f³, f_xx, f_yy


def test_true_allen_cahn_pdefind_3d():
    result = true_allen_cahn_pdefind(3, eps=0.15)
    assert result["f_xx"] == pytest.approx(0.0225)
    assert result["f_yy"] == pytest.approx(0.0225)
    assert result["f_zz"] == pytest.approx(0.0225)
    assert len(result) == 5  # f, f³, f_xx, f_yy, f_zz


def test_true_cahn_hilliard_pdefind_2d():
    result = true_cahn_hilliard_pdefind(2, eps=0.2)
    # Pure 4th
    assert result["f_xxxx"] == pytest.approx(-0.04)
    assert result["f_yyyy"] == pytest.approx(-0.04)
    # Mixed 4th
    assert result["f_xxyy"] == pytest.approx(-0.08)
    # Pure 2nd
    assert result["f_xx"] == 1.0
    assert result["f_yy"] == 1.0
    # f² * pure 2nd
    assert result["f²·f_xx"] == -3.0
    assert result["f²·f_yy"] == -3.0
    assert len(result) == 7


def test_true_allen_cahn_pdefind_raises_invalid_ndim():
    with pytest.raises(ValueError, match="ndim must be 1-8"):
        true_allen_cahn_pdefind(0)
    with pytest.raises(ValueError, match="ndim must be 1-8"):
        true_allen_cahn_pdefind(9)


# ---------------------------------------------------------------------------
# Tests — truth functions: WSINDy basis
# ---------------------------------------------------------------------------

def test_true_pme_wsindy():
    desc = _MOCK_DESCRIPTIONS_WSINDY
    w = true_pme_wsindy(desc)
    assert w[find_op(desc, '(f\u00b3)_{xx}')] == pytest.approx(1.0 / 3.0)
    assert w[find_op(desc, '(f\u00b3)_{yy}')] == pytest.approx(1.0 / 3.0)
    assert np.count_nonzero(w) == 2


# ---------------------------------------------------------------------------
# Tests — EQUATION_META
# ---------------------------------------------------------------------------

def test_equation_meta_keys():
    """EQUATION_META has the expected canonical keys."""
    expected_keys = {"ac", "sh", "ch", "pme", "ks_cell", "ks_wchaos", "ks_chaos"}
    assert set(EQUATION_META.keys()) == expected_keys


def test_equation_meta_truth_functions():
    """Each EQUATION_META entry maps to the correct truth function."""
    assert EQUATION_META["ac"][0] is true_allen_cahn
    assert EQUATION_META["sh"][0] is true_swift_hohenberg
    assert EQUATION_META["ch"][0] is true_cahn_hilliard
    assert EQUATION_META["pme"][0] is true_pme
    for key in ("ks_cell", "ks_wchaos", "ks_chaos"):
        assert EQUATION_META[key][0] is true_ks


def test_equation_meta_skip_fractions():
    assert EQUATION_META["ac"][1] == 0.0
    assert EQUATION_META["sh"][1] == 0.1
    assert EQUATION_META["ch"][1] == 0.0
    assert EQUATION_META["pme"][1] == 0.0
    for key in ("ks_cell", "ks_wchaos", "ks_chaos"):
        assert EQUATION_META[key][1] == 0.3


def test_equation_meta_masks():
    assert EQUATION_META["ac"][2] is None
    assert EQUATION_META["sh"][2] is None
    assert EQUATION_META["ch"][2] is None
    assert EQUATION_META["pme"][2] is mask_pme


def test_equation_meta_truth_outputs_match_direct():
    """Calling EQUATION_META truth fn gives same result as direct call."""
    desc = _MOCK_DESCRIPTIONS_EQUIVARIANT
    for key in ("ac", "sh", "ch", "pme"):
        truth_fn = EQUATION_META[key][0]
        direct_fn = {"ac": true_allen_cahn, "sh": true_swift_hohenberg,
                     "ch": true_cahn_hilliard, "pme": true_pme}[key]
        np.testing.assert_array_equal(truth_fn(desc), direct_fn(desc))


# ---------------------------------------------------------------------------
# Tests — normalize_eq_name
# ---------------------------------------------------------------------------

def test_normalize_eq_name_legacy():
    assert normalize_eq_name("allen_cahn") == "ac"
    assert normalize_eq_name("swift_hohenberg") == "sh"
    assert normalize_eq_name("cahn_hilliard") == "ch"
    assert normalize_eq_name("porous_medium") == "pme"
    assert normalize_eq_name("ks_cellular") == "ks_cell"
    assert normalize_eq_name("ks_weak_chaotic") == "ks_wchaos"
    assert normalize_eq_name("ks_chaotic") == "ks_chaos"


def test_normalize_eq_name_canonical():
    """Canonical names are unchanged."""
    for key in ("ac", "sh", "ch", "pme", "ks_cell", "ks_wchaos", "ks_chaos"):
        assert normalize_eq_name(key) == key


def test_normalize_eq_name_unknown_passthrough():
    assert normalize_eq_name("bogus") == "bogus"
    assert normalize_eq_name("") == ""


def test_normalize_eq_name_alias_dict():
    """_EQ_NAME_ALIASES has no dangling references — every value is a key."""
    for value in _EQ_NAME_ALIASES.values():
        assert value in _EQ_NAME_ALIASES, f"value {value!r} not a key in aliases"
