"""Shared equation configurations for sweep scripts.

Provides canonical ground-truth coefficient vectors, operator-finding
utilities, and per-equation metadata used across all sweep and experiment
scripts.  These functions are character-for-character identical across
scripts — extracting them avoids maintaining 5–6 redundant copies.
"""

import numpy as np


def find_op(descriptions, name):
    """Return column index of operator *name* in the description list.

    Parameters
    ----------
    descriptions : list of str
        Operator descriptions from the library builder.
    name : str
        Exact description string (e.g. ``'Δf'``, ``'f³'``).

    Returns
    -------
    int

    Raises
    ------
    KeyError
        If *name* is not found.
    """
    try:
        return descriptions.index(name)
    except ValueError:
        raise KeyError(f"Operator {name!r} not in library. Available: {descriptions}")


def true_allen_cahn(descriptions, eps=0.15):
    """Ground-truth coefficients for the Allen-Cahn equation.

    .. math::
        f_t = f - f^3 + \\nu\\,\\Delta f, \\qquad \\nu = \\varepsilon^2.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'f')] = 1.0
    w[find_op(descriptions, 'f³')] = -1.0
    w[find_op(descriptions, 'Δf')] = eps ** 2
    return w


def true_swift_hohenberg(descriptions, r=0.3):
    """Ground-truth coefficients for the Swift-Hohenberg equation.

    .. math::
        f_t = -(1+\\Delta)^2 f + \\alpha f - f^3,
        \\qquad \\alpha = r.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'f')] = r - 1.0
    w[find_op(descriptions, 'Δf')] = -2.0
    w[find_op(descriptions, 'Δ²f')] = -1.0
    w[find_op(descriptions, 'f³')] = -1.0
    return w


def true_ks(descriptions):
    """Ground-truth coefficients for the Kuramoto-Sivashinsky equation.

    .. math::
        f_t = -\\Delta f - \\Delta^2 f - \\|\\nabla f\\|^2.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, '‖∇f‖²')] = -0.5
    w[find_op(descriptions, 'Δf')] = -1.0
    w[find_op(descriptions, 'Δ²f')] = -1.0
    return w


def true_cahn_hilliard(descriptions, eps=0.2):
    """Ground-truth coefficients for the Cahn-Hilliard equation.

    .. math::
        f_t = -\\Delta^2 f + \\Delta(f^3 - f)
            = -\\Delta f + 3f^2\\Delta f + 6f\\|\\nabla f\\|^2
              - \\varepsilon^2 \\Delta^2 f.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'Δf')] = -1.0
    w[find_op(descriptions, 'f²·Δf')] = 3.0
    w[find_op(descriptions, 'f·‖∇f‖²')] = 6.0
    w[find_op(descriptions, 'Δ²f')] = -(eps ** 2)
    return w


def true_pme(descriptions):
    """Ground-truth coefficients for the Porous Medium Equation (m=2).

    .. math::
        f_t = \\nabla\\cdot(f^2\\nabla f) = f^2\\Delta f + 2f\\|\\nabla f\\|^2.
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, 'f²·Δf')] = 1.0
    w[find_op(descriptions, 'f·‖∇f‖²')] = 2.0
    return w


def true_pme_wsindy(descriptions):
    """Ground-truth coefficients for PME in the WSINDy basis ∂^α(f^q).

    In WSINDy, PME is exactly representable as:

        f_t = (1/3) ∂_x²(f³) + (1/3) ∂_y²(f³)

    because ∂_x²(f³) = 6 f·f_x² + 3 f²·f_{xx} and similarly for ∂_y²,
    so the sum recovers f²·Δf + 2 f·‖∇f‖².
    """
    w = np.zeros(len(descriptions))
    w[find_op(descriptions, '(f³)_{xx}')] = 1.0 / 3.0
    w[find_op(descriptions, '(f³)_{yy}')] = 1.0 / 3.0
    return w


def mask_pme(field):
    """Column mask for the Porous Medium Equation: discard near-zero values.

    Only spatial points where ``field > 0.05`` are retained.
    """
    return field.ravel() > 0.05


# Equation registry: maps name prefix → (true_fn, skip_fraction, mask_fn).
# KS variants (ks_cellular, ks_weak_chaotic, ks_chaotic) all use true_ks.
EQUATION_META = {
    "allen_cahn":       (true_allen_cahn,       0.0,  None),
    "swift_hohenberg":  (true_swift_hohenberg,  0.1,  None),
    "cahn_hilliard":    (true_cahn_hilliard,     0.0,  None),
    "porous_medium":    (true_pme,               0.0,  mask_pme),
    "ks_cellular":      (true_ks,                0.3,  None),
    "ks_weak_chaotic":  (true_ks,                0.3,  None),
    "ks_chaotic":       (true_ks,                0.3,  None),
}
