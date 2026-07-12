"""
Utility functions for equivariant PDE discovery.

Provides:
  _apply_prefilter       — periodic Gaussian prefilter (shared by spectral/polynomial backends)
  print_pde               — pretty-print a recovered PDE from a coefficient vector
  rotation_90             — rotate a 2-D periodic scalar field 90° CCW
  make_spatial_grid       — build (X, Y) coordinate grids centred at origin
  build_noise_profile     — construct a spatially-varying noise amplitude profile
  normalize_profile       — normalise a profile so mean(P²)=1 (preserves total noise power)
  apply_anisotropic_noise — apply space-dependent observation noise to snapshots
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


# ── anisotropic / space-dependent noise ─────────────────────────────────────


def make_spatial_grid(N, Lx, Ly):
    """Return (X, Y) coordinate grids centred at the origin.

    Parameters
    ----------
    N : int
        Grid resolution (square: N × N).
    Lx, Ly : float
        Domain lengths.  The grid spans [-Lx/2, Lx/2) × [-Ly/2, Ly/2).

    Returns
    -------
    X, Y : ndarray of shape (N, N)
        Meshgrid with ``indexing='ij'``.
    """
    x = np.linspace(-Lx / 2.0, Lx / 2.0, N, endpoint=False)
    y = np.linspace(-Ly / 2.0, Ly / 2.0, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    return X, Y


def build_noise_profile(X, Y, *, profile_type="gaussian_decay",
                        direction="radial", length_scale=None,
                        amplitude=1.0):
    """Build a spatially-varying noise-amplitude profile P(x, y).

    Parameters
    ----------
    X, Y : ndarray of shape (N, N)
        Coordinate grids from ``make_spatial_grid``.
    profile_type : str
        ``"gaussian_decay"`` → :math:`\\exp(-r^2/(2s^2))`  (noise maximal at
        origin, decays outward).  ``"gaussian_growth"`` →
        :math:`1 + A\\,(1 - \\exp(-r^2/(2s^2)))`  (baseline noise at origin,
        grows to :math:`1+A` at large :math:`r`).
    direction : str
        ``"radial"`` → :math:`r = \\sqrt{x^2 + y^2}`  (O(2)-symmetric).
        ``"axial"`` → :math:`r = |x|`  (breaks rotational symmetry).
    length_scale : float or None
        Gaussian width :math:`s` in physical units.  Default: ``Lx / 4``
        (inferred from X range).
    amplitude : float
        Growth amplitude :math:`A` (only used for ``profile_type="gaussian_growth"``).
        Default 1.0 → noise doubles at large :math:`r`.

    Returns
    -------
    P : ndarray of shape (N, N)
        Raw (unnormalised) profile.
    """
    if direction == "radial":
        R = np.sqrt(X ** 2 + Y ** 2)
    elif direction == "axial":
        R = np.abs(X)
    else:
        raise ValueError(f"Unknown direction: {direction!r}")

    if length_scale is None:
        length_scale = _infer_length_scale(X)

    if profile_type == "gaussian_decay":
        P = np.exp(-R ** 2 / (2.0 * length_scale ** 2))
    elif profile_type == "gaussian_growth":
        P = 1.0 + amplitude * (1.0 - np.exp(-R ** 2 / (2.0 * length_scale ** 2)))
    else:
        raise ValueError(f"Unknown profile_type: {profile_type!r}")

    return P


def _infer_length_scale(X):
    """Infer Lx from the X coordinate grid."""
    x_vals = X[:, 0]
    dx = x_vals[1] - x_vals[0]
    Lx = x_vals[-1] - x_vals[0] + dx
    return Lx / 4.0


def normalize_profile(P):
    """Normalise a noise profile so that ``mean(P²) = 1``.

    This preserves the **total noise power** relative to the isotropic
    (homogeneous) case: for the same ``obs_noise`` fraction, the expected
    Frobenius norm of the noise field is unchanged.

    Parameters
    ----------
    P : ndarray of shape (N, N)
        Raw profile from ``build_noise_profile``.

    Returns
    -------
    P_norm : ndarray of shape (N, N)
        Normalised profile with ``mean(P_norm ** 2) ≈ 1``.
    """
    rms = np.sqrt(np.mean(P ** 2))
    return P / rms


def apply_anisotropic_noise(snapshots, obs_noise, profile, rng):
    """Apply space-dependent observation noise to a stack of clean snapshots.

    The noise at each grid point :math:`(i, j)` and snapshot :math:`t` is

    .. math::

        \\sigma \\cdot P_{ij} \\cdot \\mathcal{N}(0, 1)

    where :math:`\\sigma = \\text{obs\\_noise} \\cdot \\text{std}(\\text{snapshots})`
    (computed over all snapshots and all spatial points, identical to the
    isotropic convention).

    Parameters
    ----------
    snapshots : ndarray of shape (T, N, N)
        Clean field snapshots.
    obs_noise : float
        Noise fraction (e.g. 0.02 for 2 %).
    profile : ndarray of shape (N, N)
        Normalised noise-amplitude profile (e.g. from ``normalize_profile``).
    rng : numpy.random.Generator
        Seeded random generator for reproducibility.

    Returns
    -------
    noisy : ndarray of shape (T, N, N)
        Snapshots with anisotropic observation noise added.
    """
    sigma = obs_noise * snapshots.std()
    noise = sigma * profile[np.newaxis, :, :] * rng.standard_normal(snapshots.shape)
    return snapshots + noise


# ── spectral low-pass filter ────────────────────────────────────────────────


def spectral_lowpass(f, Lx, Ly, cutoff):
    """Apply a 2-D isotropic spectral low-pass filter.

    Zeroes all Fourier modes with physical wavenumber magnitude |k| > cutoff
    and returns the result in space domain.  The filter is O(2)-symmetric
    (isotropic cutoff sphere), so it commutes with all O(2)-equivariant
    operators in the library.

    Parameters
    ----------
    f      : ndarray, shape (Nx, Ny)
    Lx, Ly : float - domain lengths
    cutoff : float - maximum physical wavenumber to retain

    Returns
    -------
    f_filtered : ndarray, shape (Nx, Ny), real
    """
    Nx, Ny = f.shape
    kx = 2.0 * np.pi / Lx * np.fft.fftfreq(Nx) * Nx
    ky = 2.0 * np.pi / Ly * np.fft.fftfreq(Ny) * Ny
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    fh = np.fft.fft2(f)
    fh[KX ** 2 + KY ** 2 > cutoff ** 2] = 0.0
    return np.real(np.fft.ifft2(fh))
