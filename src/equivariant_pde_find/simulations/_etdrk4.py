import numpy as np


def _etdrk4_coeffs(Lk, dt, M=32, r=1.0):
    """Precompute ETDRK4 coefficients for a diagonal linear operator Lk.

    Uses the Kassam-Trefethen (2005) contour-integral method to evaluate
    the phi-functions stably near z=0 (avoids catastrophic cancellation for
    small or zero wavenumbers).

    Parameters
    ----------
    Lk : ndarray
        Diagonal of the linear operator in Fourier space (any shape).
    dt : float
        Time step.
    M : int
        Number of contour quadrature points (default 32; 16 is also fine).
    r : float
        Contour radius (default 1.0).

    Returns
    -------
    E, E2, phi1_dt2, f1, f2, f4 : ndarrays of same shape as Lk
        E        = exp(Lk * dt)
        E2       = exp(Lk * dt/2)
        phi1_dt2 = (dt/2) * phi1(Lk * dt/2)  -- coefficient for stages
        f1, f2, f4               -- coefficients for the final update
    """
    # Contour points, shape (M,)
    contour = r * np.exp(2j * np.pi * np.arange(M) / M)

    # Half- and full-step z, shape (*Lk.shape, 1) for broadcasting
    z2 = Lk[..., np.newaxis] * (dt / 2)       # Lk * dt/2 + (M,)
    z  = Lk[..., np.newaxis] * dt             # Lk * dt   + (M,)

    zeta2 = z2 + contour  # shape (*Lk.shape, M)
    zeta  = z  + contour

    E2 = np.exp(Lk * (dt / 2))
    E  = np.exp(Lk * dt)

    # phi1 at half step:  (dt/2) * (e^{z2} - 1) / z2
    # Limiting value as z2 -> 0 is 1.  Handle the (rare) case where a
    # quadrature point lands exactly on 0 to avoid 0/0 NaN.
    _safe_z2 = np.where(np.abs(zeta2) < 1e-10, 1.0, zeta2)
    _phi1_integrand = np.where(
        np.abs(zeta2) < 1e-10, 1.0, (np.exp(zeta2) - 1) / _safe_z2
    )
    phi1_dt2 = (dt / 2) * np.mean(_phi1_integrand.real, axis=-1)

    # Final-step coefficients (Cox-Matthews formulae, Kassam-Trefethen Eq. 2.4)
    # All three functions have the same limit dt/6 as z -> 0.
    ez = np.exp(zeta)
    _safe_z3 = np.where(np.abs(zeta) < 1e-10, 1.0, zeta ** 3)
    _lim = dt / 6  # common limit of all three integrands * dt as z->0

    def _coeff(numerator):
        integrand = np.where(np.abs(zeta) < 1e-10, _lim, dt * numerator / _safe_z3)
        return np.mean(integrand.real, axis=-1)

    f1 = _coeff(-4 - zeta + ez * (4 - 3 * zeta + zeta ** 2))
    f2 = _coeff(2 + zeta + ez * (-2 + zeta))
    f4 = _coeff(-4 - 3 * zeta - zeta ** 2 + ez * (4 - zeta))

    return E, E2, phi1_dt2, f1, f2, f4


def _etdrk4_run(f0, nonlin_hat_fn, E, E2, phi1_dt2, f1, f2, f4, n_steps, save_every):
    """ETDRK4 time loop.  Works entirely in Fourier space.

    Parameters
    ----------
    f0 : ndarray (physical space, shape (N, N))
    nonlin_hat_fn : callable
        Maps Fourier-space array û -> FFT of the nonlinear term N(u).
    E, E2, phi1_dt2, f1, f2, f4 : precomputed coefficient arrays
    n_steps, save_every : int

    Returns
    -------
    snapshots : ndarray of shape (n_saved, N, N)  [physical space]
    times : ndarray of shape (n_saved,)
    """
    uhat = np.fft.fft2(f0)
    snapshots, steps = [], []

    for step in range(n_steps):
        Nu  = nonlin_hat_fn(uhat)
        a   = E2 * uhat + phi1_dt2 * Nu
        Na  = nonlin_hat_fn(a)
        b   = E2 * uhat + phi1_dt2 * Na
        Nb  = nonlin_hat_fn(b)
        c   = E2 * a   + phi1_dt2 * (2 * Nb - Nu)
        Nc  = nonlin_hat_fn(c)
        uhat = E * uhat + f1 * Nu + 2 * f2 * (Na + Nb) + f4 * Nc

        if step % save_every == 0:
            snapshots.append(np.real(np.fft.ifft2(uhat)))
            steps.append(step)

    return np.array(snapshots), np.array(steps, dtype=float)
