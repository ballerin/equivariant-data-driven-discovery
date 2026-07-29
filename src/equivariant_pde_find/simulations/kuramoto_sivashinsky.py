"""Kuramoto-Sivashinsky equation simulator (2-D only; see the ``ndim``
parameter's docstring for why the higher-dimensional generalization done
for the other equations does not extend here without further validation).
"""

import numpy as np

from ._imex import _run
from ._etdrk4 import _etdrk4_coeffs, _etdrk4_run


def kuramoto_sivashinsky(
    N=64,
    L=16 * np.pi,
    T_end=400.0,
    dt=0.25,
    save_every=1,
    f0=None,
    seed=42,
    obs_noise=0.0,
    dyn_noise=0.0,
    scheme="etdrk4",
    dealias=True,
    ndim=2,
):
    """Simulate the 2D Kuramoto–Sivashinsky equation via pseudo-spectral methods.

    PDE:
        f_t = -1/2 |∇f|^2 - Δf - Δ^2 f

    The linear Fourier symbol is L_k = K^2 - K^4.

    Parameters
    ----------
    N : int
        Grid points per direction (power of 2 recommended).
    L : float
        Domain side length (periodic, [0, L]^2).
    T_end : float
        End time.
    dt : float
        Time step.  Default 0.25 (suitable for ETDRK4; use 0.05 for IMEX).
    save_every : int
        Save a snapshot every this many steps.
    f0 : ndarray of shape (N, N) or None
        Initial condition.
    seed : int
        RNG seed.
    obs_noise : float
        Observation noise as a fraction of ``clean_snapshots.std()``.
    dyn_noise : float
        Dynamical noise (only available with ``scheme='imex'``).
    scheme : {'etdrk4', 'imex'}
        Time-stepping scheme.  'etdrk4' (default) uses the fourth-order
        exponential Runge-Kutta method of Cox & Matthews (2002) with the
        contour-integral fix of Kassam & Trefethen (2005).  'imex' uses the
        first-order implicit-explicit Euler scheme (legacy, kept for
        reproducibility).
    dealias : bool
        Apply the 2/3-rule dealiasing mask before computing the nonlinear
        term.  Strongly recommended (default True); avoids quadratic
        aliasing that contaminates the regression library.
    ndim : int
        Number of spatial dimensions.  Only ``ndim=2`` is currently
        supported: the ETDRK4 coefficients and time-stepping loop here
        are hardcoded around ``fft2``/``ifft2`` and have not been
        verified for other dimensions (the exponential-integrator
        stability analysis is more delicate than the plain IMEX schemes
        used elsewhere).  Any other value raises ``ValueError``.

    Returns
    -------
    snapshots : ndarray of shape (n_saved, N, N)
    times : ndarray of shape (n_saved,)
    """
    if ndim != 2:
        raise ValueError(
            f"kuramoto_sivashinsky only supports ndim=2 (got ndim={ndim}); "
            "the ETDRK4/IMEX integrators here are hardcoded for 2-D FFTs "
            "and have not been verified for other dimensions."
        )

    rng = np.random.default_rng(seed)

    k = np.fft.fftfreq(N) * N * (2 * np.pi / L)
    KX, KY = np.meshgrid(k, k, indexing="ij")
    K2 = KX**2 + KY**2
    K4 = K2**2

    Lk = K2 - K4  # linear Fourier symbol

    # 2/3-rule dealiasing mask (zero top 1/3 of modes in each direction)
    if dealias:
        k_max = (N // 2) * (2 * np.pi / L)
        dealias_mask = (np.abs(KX) <= 2 / 3 * k_max) & (np.abs(KY) <= 2 / 3 * k_max)
    else:
        dealias_mask = None

    if f0 is None:
        f0 = 0.1 * rng.standard_normal((N, N))

    n_steps = int(round(T_end / dt)) + 1

    if scheme == "etdrk4":
        E, E2, phi1_dt2, f1, f2, f4 = _etdrk4_coeffs(Lk, dt)

        def nonlin_hat(uhat):
            """FFT of the nonlinear term -1/2 |∇u|^2, given û in Fourier space."""
            if dealias_mask is not None:
                uhat = uhat * dealias_mask
            fx = np.real(np.fft.ifft2(1j * KX * uhat))
            fy = np.real(np.fft.ifft2(1j * KY * uhat))
            return np.fft.fft2(-0.5 * (fx * fx + fy * fy))

        snapshots, steps = _etdrk4_run(
            f0, nonlin_hat, E, E2, phi1_dt2, f1, f2, f4, n_steps, save_every
        )
        times = steps * dt

    else:  # imex
        denom = 1.0 - dt * Lk

        def nonlin(f):
            fhat = np.fft.fft2(f)
            if dealias_mask is not None:
                fhat = fhat * dealias_mask
            fx = np.real(np.fft.ifft2(1j * KX * fhat))
            fy = np.real(np.fft.ifft2(1j * KY * fhat))
            return -0.5 * (fx * fx + fy * fy)

        snapshots, times = _run(f0, nonlin, denom, dt, n_steps, save_every, dyn_noise, rng)

    if obs_noise > 0.0:
        sigma = obs_noise * snapshots.std()
        snapshots = snapshots + sigma * rng.standard_normal(snapshots.shape)

    return snapshots, times
