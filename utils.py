"""
Utility functions for equivariant PDE discovery.

Provides:
  print_pde    — pretty-print a recovered PDE from a coefficient vector
  rotation_90  — rotate a 2-D periodic scalar field 90° CCW
"""

import numpy as np


def print_pde(w, descriptions, lhs="f_t"):
    """Pretty-print an identified PDE from coefficient vector w."""
    terms = [(c, d) for c, d in zip(w, descriptions) if abs(c) > 1e-10]
    if not terms:
        print(f"{lhs} = 0")
        return
    rhs = " + ".join(f"({c:+.5f}){d}" for c, d in terms)
    print(f"{lhs} = {rhs}")


def rotation_90(f):
    """Rotate field 90° CCW on a square periodic grid.

    For f on an (N, N) grid: R[f][i, j] = f[j, N-1-i].
    """
    return np.rot90(f, k=1)
