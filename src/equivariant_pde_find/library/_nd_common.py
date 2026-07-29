"""
Shared N-D helpers for the coordinate-basis library builders
(pdefind.py, polyjet.py, wsindy.py): multi-index enumeration and
per-field derivative computation.

Not part of the public API -- internal to the library/ subpackage.
"""

from ..derivatives.polynomial import poly_partial_nd


def indices_of_order(ndim, total):
    """All ndim-tuples of nonneg ints summing to exactly `total`, in the
    historical 'leading axis varies slowest' lexicographic order."""
    if ndim == 1:
        return [(total,)]
    out = []
    for a in range(total + 1):
        for rest in indices_of_order(ndim - 1, total - a):
            out.append((a,) + rest)
    return out


def enumerate_multi_indices(ndim, max_order):
    """All ndim-tuples of nonneg ints with sum <= max_order, grouped by
    total order ascending (0, 1, 2, ...) -- matches the historical 2-D/3-D
    (a, b[, c]) enumeration order exactly."""
    result = []
    for total in range(max_order + 1):
        result.extend(indices_of_order(ndim, total))
    return result


def compute_derivs_nd(f, Ls, multi_indices, method, zero_value="f"):
    """Compute {multi_index: ndarray} for f via the given backend method
    ('spectral' or 'polynomial'), sharing one forward FFT for spectral.

    Parameters
    ----------
    zero_value : 'f' or 'ones'
        What the zero multi-index (the "0th derivative") maps to.
        'f' -- the field itself (jet-variable convention: f is jet
        variable of order 0; used by polyjet.py/wsindy.py).
        'ones' -- the constant 1 (PDE-FIND convention, where the f^p
        polynomial prefactor is tracked separately from the derivative
        factor; used by pdefind.py).
    """
    import numpy as np

    nonzero_indices = [mi for mi in multi_indices if any(mi)]
    if method == "spectral":
        from ..derivatives.spectral import spectral_partials_nd_batch
        out = spectral_partials_nd_batch(f, Ls, nonzero_indices) if nonzero_indices else {}
    else:
        out = {mi: poly_partial_nd(f, Ls, mi) for mi in nonzero_indices}

    zero_mi = tuple(0 for _ in range(f.ndim))
    if zero_mi in multi_indices:
        out[zero_mi] = f.copy() if zero_value == "f" else np.ones_like(f)
    return out
