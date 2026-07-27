import numpy as np


def make_regression_data(snapshots, times, n_per_snap=200, skip_fraction=0.0, seed=0):
    """
    Extract (f_snaps, ft_snaps, sample_indices) from a snapshot array.

    Uses central finite differences in time to estimate f_t.  Optionally
    skips an initial transient fraction of the trajectory.

    Parameters
    ----------
    snapshots : ndarray of shape (n_saved, N, N)
    times : ndarray of shape (n_saved,)
    n_per_snap : int
        Number of spatial points sampled per snapshot.
    skip_fraction : float in [0, 1)
        Skip this fraction of snapshots from the start (transient removal).
    seed : int
        RNG seed for spatial subsampling.

    Returns
    -------
    f_snaps : ndarray of shape (n_interior, N, N)
        Field snapshots at interior time points.
    ft_snaps : ndarray of shape (n_interior, N, N)
        Time derivatives estimated by central finite difference.
    idx_sets : list of ndarray
        Spatial indices sampled for each snapshot.
    dt_snap : float
        Time between consecutive saved snapshots.
    """
    skip = max(0, int(skip_fraction * len(snapshots)))
    snapshots = snapshots[skip:]
    times = times[skip:]

    dt_snap = times[1] - times[0]
    ft_snaps = (snapshots[2:] - snapshots[:-2]) / (2.0 * dt_snap)
    f_snaps = snapshots[1:-1]

    N_total = np.prod(f_snaps.shape[1:]).item()
    rng = np.random.default_rng(seed)
    idx_sets = [rng.choice(N_total, n_per_snap, replace=False) for _ in range(len(f_snaps))]

    return f_snaps, ft_snaps, idx_sets, dt_snap
