"""
Change-of-basis matrices between the original multigraph basis and the
modified multigraph (IBP-friendly).
"""

import numpy as np

from .multigraph_operators import (
    get_mod_multigraph_operators,
    get_operators,
)


# -----------------------------------------------------------------------
# Block matrices mapping mod. multigraph → original (w_orig = M @ w_mod)
# -----------------------------------------------------------------------

# Each block corresponds to a specific (p, n) subspace.  For subspaces
# where the two bases share the same operators, the block is identity
# The following matrices encode the mathematical identities that relate
# the two bases, as described in the reference paper.

# (p=1, n=2):
# - Δf² → f·Δf
# - shared ‖∇f‖²
_M_12 = np.array([[2.0, 0.0],
                  [2.0, 1.0]])

# (p=1, n=3): 
# - Δf³ → f²·Δf  
# - shared f·‖∇f‖²
_M_13 = np.array([[3.0, 0.0],
                  [6.0, 1.0]])

# (p=2, n=2):
# - Δ²f² → f·Δ²f
# - Δ‖∇f‖² → ⟨∇Δf,∇f⟩
# - shared ‖∇²f‖²
# - shared (Δf)²
_M_22 = np.array([[2.0, 0.0, 0.0, 0.0],
                  [8.0, 2.0, 0.0, 0.0],
                  [4.0, 2.0, 1.0, 0.0],
                  [2.0, 0.0, 0.0, 1.0]])

# (p=2, n=3):
# Δ²f³→f²·Δ²f
# Δ(f·‖∇f‖²)→f·⟨∇Δf,∇f⟩
# - shared f·‖∇²f‖²
# - shared f·Δf·Δf
# - shared ∇²f(∇f,∇f)
# - shared ‖∇f‖²·Δf
_M_23 = np.array([[3.0,  0.0, 0.0, 0.0, 0.0, 0.0],
                  [24.0, 2.0, 0.0, 0.0, 0.0, 0.0],
                  [12.0, 2.0, 1.0, 0.0, 0.0, 0.0],
                  [6.0,  0.0, 0.0, 1.0, 0.0, 0.0],
                  [24.0, 4.0, 0.0, 0.0, 1.0, 0.0],
                  [12.0, 1.0, 0.0, 0.0, 0.0, 1.0]])


# -----------------------------------------------------------------------
# Index lookup helpers
# -----------------------------------------------------------------------

def _build_index_map(descs_mod, descs_orig):
    """Build a dict mapping operator description → column index."""
    return {desc: i for i, desc in enumerate(descs_mod)}, \
           {desc: i for i, desc in enumerate(descs_orig)}


def _get_block_indices(descs_mod, descs_orig, mod_names, orig_names):
    """Return (idx_mod, idx_orig) index slices for a (p,n) block."""
    mod_map, orig_map = _build_index_map(descs_mod, descs_orig)
    idx_mod = [mod_map[n] for n in mod_names]
    idx_orig = [orig_map[n] for n in orig_names]
    return idx_mod, idx_orig


# -----------------------------------------------------------------------
# Block specifications: for each (p, n) that differs between bases,
# list the string identifiers of the differential operators in the SAME
# order as the block matrices above.
#
# Triple of (mod_names, orig_names, block_matrix)
# - mod_names: list of operator descriptions in the mod. multigraph basis
# - orig_names: list of operator descriptions in the og. multigraph basis
# - block_matrix: the change-of-basis matrix for this (p, n) subspace
# -----------------------------------------------------------------------

_BLOCKS = [
    # (p=1, n=2):
    # - Δf² → f·Δf
    # - shared ‖∇f‖²
    (["Δf²", "‖∇f‖²"],
     ["f·Δf", "‖∇f‖²"],
     _M_12),

    # (p=1, n=3): 
    # - Δf³ → f²·Δf  
    # - shared f·‖∇f‖²
    (["Δf³", "f·‖∇f‖²"],
     ["f²·Δf", "f·‖∇f‖²"],
     _M_13),

    # (p=2, n=2):
    # - Δ²f²→f·Δ²f
    # - Δ‖∇f‖²→⟨∇Δf,∇f⟩
    # - shared ‖∇²f‖²
    # - shared Δf·Δf
    (["Δ²f²", "Δ‖∇f‖²", "‖∇²f‖²", "Δf·Δf"],
     ["f·Δ²f", "⟨∇Δf,∇f⟩", "‖∇²f‖²", "Δf·Δf"],
     _M_22),

    # (p=2, n=3):
    # Δ²f³→f²·Δ²f
    # Δ(f·‖∇f‖²)→f·⟨∇Δf,∇f⟩
    # - shared f·‖∇²f‖²
    # - shared f·Δf·Δf
    # - shared ∇²f(∇f,∇f)
    # - shared ‖∇f‖²·Δf
    (["Δ²f³", "Δ(f·‖∇f‖²)", "f·‖∇²f‖²", "f·Δf·Δf",
      "∇²f(∇f,∇f)", "‖∇f‖²·Δf"],
     ["f²·Δ²f", "f·⟨∇Δf,∇f⟩", "f·‖∇²f‖²", "f·Δf·Δf",
      "∇²f(∇f,∇f)", "‖∇f‖²·Δf"],
     _M_23),
]


# -----------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------


def mod_multigraph_to_multigraph_matrix(max_p=2, max_n=3):
    """Build the change-of-basis matrix M such that w_orig = M @ w_mod.

    For most (p, n) subspaces the two bases are identical, so the
    corresponding diagonal block of M is an identity matrix.  For the
    subspaces where the modified basis is genuinely different, the blocks
    are filled with the coefficients from the reference paper.

    Parameters
    ----------
    max_p, max_n : int - passed to get_operators / get_mod_multigraph_operators

    Returns
    -------
    M : ndarray (K, K) where K = len(get_operators(max_p, max_n))
    """
    descs_orig = [d for _, d in get_operators(max_p, max_n)]
    descs_mod = [d for _, d in get_mod_multigraph_operators(max_p, max_n)]
    K = len(descs_orig)
    M = np.eye(K)

    for mod_names, orig_names, block in _BLOCKS:
        # Build index slices for this block
        mod_map = {d: i for i, d in enumerate(descs_mod)}
        orig_map = {d: i for i, d in enumerate(descs_orig)}
        try:
            idx_mod = [mod_map[n] for n in mod_names]
            idx_orig = [orig_map[n] for n in orig_names]
        except KeyError as e:
            # Block not present in this (max_p, max_n) configuration
            continue
        for i_local, i_full in enumerate(idx_orig):
            for j_local, j_full in enumerate(idx_mod):
                M[i_full, j_full] = block[i_local, j_local]

    return M


def multigraph_to_mod_multigraph_matrix(max_p=2, max_n=3):
    """Build M⁻¹ s.t. w_mod = M⁻¹ @ w_orig.
    This is simply the inverse of mod_multigraph_to_multigraph_matrix(max_p, max_n).
    """
    M = mod_multigraph_to_multigraph_matrix(max_p, max_n)
    return np.linalg.inv(M)


def convert_coefficients(w, from_basis, to_basis, max_p=2, max_n=3):
    """Convert a coefficient vector between bases.

    Parameters
    ----------
    w : ndarray (K,)
        Coefficient vector in ``from_basis``.
    from_basis : str
        One of ``'multigraph'`` or ``'mod_multigraph'``.
    to_basis : str
        One of ``'multigraph'`` or ``'mod_multigraph'``.
    max_p, max_n : int

    Returns
    -------
    w_out : ndarray (K,)
        Coefficient vector in ``to_basis``.
    """
    if from_basis == to_basis:
        return w.copy()
    if from_basis == "mod_multigraph" and to_basis == "multigraph":
        M = mod_multigraph_to_multigraph_matrix(max_p, max_n)
        return M @ w
    if from_basis == "multigraph" and to_basis == "mod_multigraph":
        M_inv = multigraph_to_mod_multigraph_matrix(max_p, max_n)
        return M_inv @ w
    raise ValueError(
        f"Unknown basis conversion: {from_basis} → {to_basis}. "
        f"Expected 'multigraph' or 'mod_multigraph'."
    )
