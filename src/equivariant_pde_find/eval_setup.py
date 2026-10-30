"""Ground truth and the evaluation design matrix for one dataset.

Scoring a discovered equation is backend- and formulation-independent: every
library's coefficients are mapped analytically into a common polynomial jet
basis and compared there, against a pointwise spectral design matrix built on
*clean* data.  Setting that up takes ~90 lines and used to be inlined in
``scripts/sweep.py``'s main loop, which meant a second driver could not score
anything without copying it.  :func:`prepare_evaluation` is that block,
extracted verbatim, so ``sweep.py`` and the train/validation/test driver share
one implementation.

The train and validation tiers of a split protocol never call this: they only
need a held-out residual, which reads no ground truth.  It is the test tier
that scores against truth, so the cost is paid once per test trajectory
rather than once per dataset.
"""

import time

import numpy as np

from .equations import CH_EPSILON
from .evaluation import FactoredEval, build_eval_theta_jet
from .library.conversions import (
    _MULTIGRAPH_TO_JET_DECOMPOSITION,
    build_wsindy_to_jet_map,
    collect_needed_jet_monomials,
)
from .library import (
    build_multigraph_library,
    build_pdefind_library,
    build_wsindy_library,
)
from .simulations import make_regression_data


class Evaluation:
    """Everything needed to score a coefficient vector for one dataset.

    Attributes
    ----------
    desc_eq_ref, true_eq : the multigraph column names and the ground-truth
        coefficients in that basis.
    desc_jet, true_jet : the common jet basis and the ground truth mapped into
        it.
    M_ws_to_jet : WSINDy-to-jet map.
    Theta_eval, Ft_eval : the clean pointwise spectral system, or a
        :class:`~equivariant_pde_find.evaluation.FactoredEval` with
        ``Ft_eval=None`` when ``eval_mode="qr"``.
    f_snaps, ft_snaps, idx_sets : per-snapshot regression data (pointwise
        formulations only).
    n_snaps : number of interior snapshots.
    """

    __slots__ = ("desc_eq_ref", "true_eq", "desc_jet", "true_jet",
                 "M_ws_to_jet", "Theta_eval", "Ft_eval",
                 "f_snaps", "ft_snaps", "idx_sets", "n_snaps")

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def ground_truth(snapshot0, Lx, Ly, eq_name, true_fn, ch_eps=None):
    """Ground-truth coefficients in the multigraph basis.

    Returns ``(desc_eq_ref, true_eq)``.  The per-equation parameters match
    the values the simulators were run with.
    """
    _, desc_eq_ref = build_multigraph_library(snapshot0, Lx, Ly)
    if eq_name == "ch":
        true_eq = true_fn(desc_eq_ref,
                          eps=CH_EPSILON if ch_eps is None else ch_eps)
    elif eq_name == "ac":
        true_eq = true_fn(desc_eq_ref, eps=0.15)
    elif eq_name == "sh":
        true_eq = true_fn(desc_eq_ref, r=0.3)
    else:
        true_eq = true_fn(desc_eq_ref)
    return desc_eq_ref, true_eq


def prepare_evaluation(snapshots_proc, snapshots_clean_eval, times, Lx, Ly,
                       eq_name, true_fn, skip_frac, mask_fn=None,
                       ch_eps=None, eval_mode="dense", verbose=True):
    """Ground truth and the clean evaluation system for one dataset.

    Parameters mirror what ``sweep.py`` has in scope at this point:
    ``snapshots_proc`` is the (possibly noisy, possibly filtered) data that
    libraries are built from, ``snapshots_clean_eval`` the clean counterpart
    that scoring uses, and ``mask_fn`` an optional per-snapshot validity mask
    (PME's positivity support).

    ``eval_mode="qr"`` replaces the dense evaluation matrix by the R factor of
    its QR decomposition: ``E_coeff`` agrees to round-off and the other
    metrics are identical, but scoring no longer scales with the grid.
    """
    desc_eq_ref, true_eq = ground_truth(snapshots_proc[0], Lx, Ly, eq_name,
                                        true_fn, ch_eps)

    # Prepare regression data
    f_snaps_raw, ft_snaps_raw, idx_sets_raw, _ = make_regression_data(
        snapshots_proc, times, n_per_snap=200, skip_fraction=skip_frac, seed=0,
    )

    # Apply mask if needed
    if mask_fn is not None:
        rng = np.random.default_rng(1)
        new_idx = []
        for s in range(len(f_snaps_raw)):
            valid = np.where(mask_fn(f_snaps_raw[s]))[0]
            if len(valid) <= 200:
                new_idx.append(valid)
            else:
                new_idx.append(np.sort(
                    rng.choice(valid, size=200, replace=False)))
        f_snaps, ft_snaps, idx_sets = f_snaps_raw, ft_snaps_raw, new_idx
    else:
        f_snaps, ft_snaps, idx_sets = f_snaps_raw, ft_snaps_raw, idx_sets_raw

    n_snaps = len(f_snaps)
    if verbose:
        print(f"    interior snaps: {n_snaps}  points: {n_snaps * 200}",
              flush=True)

    # ── Jet descriptions & ground truth (backend-independent) ─
    # Build pdefind descriptions upfront (backend-independent: column names
    # are the same regardless of spectral/polynomial)
    _, desc_pdefind = build_pdefind_library(
        snapshots_proc[0], Lx, Ly, max_order=4, max_poly=3,
    )
    _, desc_ws = build_wsindy_library(
        snapshots_proc[0], Lx, Ly, max_order=4, max_poly=3,
    )
    # Dynamic jet: union of all library decompositions
    desc_jet = collect_needed_jet_monomials(
        desc_mg=desc_eq_ref, desc_pdefind=desc_pdefind, desc_ws=desc_ws,
    )
    # Ground truth in jet basis (analytical — no backend)
    true_jet = np.zeros(len(desc_jet))
    for name, coeff in zip(desc_eq_ref, true_eq):
        if abs(coeff) < 1e-12:
            continue
        if name not in _MULTIGRAPH_TO_JET_DECOMPOSITION:
            continue
        for jet_name, mult in _MULTIGRAPH_TO_JET_DECOMPOSITION[name]:
            true_jet[desc_jet.index(jet_name)] += coeff * mult

    M_ws_to_jet = build_wsindy_to_jet_map(desc_ws, desc_jet)

    # ── Build evaluation Theta (pointwise, spectral) ────
    # Always use pointwise formulation + spectral derivatives on clean
    # snapshots.  The polynomial jet basis is the universal comparison space
    # into which all library coefficients are analytically mapped.
    t_eval = time.time()
    Theta_eval, Ft_eval = build_eval_theta_jet(
        snapshots_clean_eval, times, Lx, Ly,
        skip_fraction=skip_frac, seed=0, backend="spectral",
        jet_monomials=desc_jet,
    )
    if verbose:
        print(f"    eval Theta: {Theta_eval.shape}  active_true="
              f"{int(np.sum(np.abs(true_jet) > 1e-10))}"
              f"  ({time.time() - t_eval:.1f}s)", flush=True)
    if eval_mode == "qr":
        t_eval = time.time()
        Theta_eval = FactoredEval(Theta_eval, Ft_eval)
        Ft_eval = None
        if verbose:
            print(f"    eval QR factor: {Theta_eval.R.shape}  "
                  f"({time.time() - t_eval:.1f}s)", flush=True)

    return Evaluation(
        desc_eq_ref=desc_eq_ref, true_eq=true_eq,
        desc_jet=desc_jet, true_jet=true_jet, M_ws_to_jet=M_ws_to_jet,
        Theta_eval=Theta_eval, Ft_eval=Ft_eval,
        f_snaps=f_snaps, ft_snaps=ft_snaps, idx_sets=idx_sets,
        n_snaps=n_snaps,
    )
