import numpy as np


def _imex_step(f, nonlin_fn, denom, dt, dyn_noise, rng):
    """Single IMEX step: implicit linear + explicit nonlinear + optional noise."""
    fhat = np.fft.fft2(f) + dt * np.fft.fft2(nonlin_fn(f))
    if dyn_noise > 0.0:
        xi = rng.standard_normal(f.shape)
        fhat += np.sqrt(dt) * dyn_noise * np.fft.fft2(xi)
    return np.real(np.fft.ifft2(fhat / denom))


def _run(f0, nonlin_fn, denom, dt, n_steps, save_every, dyn_noise, rng):
    """Run the time loop, returning (snapshots, step_indices)."""
    f = f0.copy()
    snapshots, steps = [], []
    for step in range(n_steps):
        f = _imex_step(f, nonlin_fn, denom, dt, dyn_noise, rng)
        if step % save_every == 0:
            snapshots.append(f.copy())
            steps.append(step)
    return np.array(snapshots), np.array(steps, dtype=float) * dt
