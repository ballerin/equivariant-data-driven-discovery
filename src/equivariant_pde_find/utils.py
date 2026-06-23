"""
Utility functions for equivariant PDE discovery.

Provides:
  _apply_prefilter — periodic Gaussian prefilter (shared by spectral/polynomial backends)
  print_pde       — pretty-print a recovered PDE from a coefficient vector
  rotation_90     — rotate a 2-D periodic scalar field 90° CCW
"""

import numpy as np


def _apply_prefilter(f, Lx, Ly, prefilter_sigma):
    """Apply a periodic Gaussian prefilter to 2-D field `f`.

    Parameters
    ----------
    f : ndarray of shape (Nx, Ny)
    Lx, Ly : float — domain lengths
    prefilter_sigma : float, (float, float), or None
        Gaussian std in physical units. If None or 0, returns f unchanged.
    """
    if not prefilter_sigma:
        return f
    Nx, Ny = f.shape
    kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
    ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    if np.isscalar(prefilter_sigma):
        sx = float(prefilter_sigma)
        sy = float(prefilter_sigma)
    else:
        sx = float(prefilter_sigma[0])
        sy = float(prefilter_sigma[1])
    fh = np.fft.fft2(f)
    G = np.exp(-0.5 * ((sx * KX) ** 2 + (sy * KY) ** 2))
    return np.real(np.fft.ifft2(fh * G))


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
