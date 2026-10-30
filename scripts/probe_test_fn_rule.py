#!/usr/bin/env python3
"""Report the Messenger & Bortz test-function parameters for a data directory.

No regression is run: this is the cheap pre-flight for a space-time weak
sweep.  It prints, per (equation, noise), the rule's (m, p, m_t, p_t), the
critical wavenumbers it derived them from, whether either support sits at its
admissibility cap, and the resulting number of space-time query points -- the
row count the sweep will have to fit.

Usage
-----
    python scripts/probe_test_fn_rule.py --data-dir data/paper
    python scripts/probe_test_fn_rule.py --data-dir data/large --tau-hat 3
"""

import argparse
import os
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(SCRIPT_DIR, "..", "src")))
sys.path.insert(0, SCRIPT_DIR)

import sweep as S                                              # noqa: E402
import equivariant_pde_find.library.weak_form as wf            # noqa: E402

NOISE_LABELS = [("clean", 0.0), ("noise2pct", 0.02), ("noise5pct", 0.05)]
EQUATIONS = ["ac", "sh", "ch", "pme", "ks_cell", "ks_wchaos", "ks_chaos"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="data/paper")
    ap.add_argument("--equations", nargs="+", default=EQUATIONS)
    ap.add_argument("--tau", type=float, default=S.TAU)
    ap.add_argument("--tau-hat", type=float, default=2.0)
    ap.add_argument("--space-stride", type=int, default=S.SPACE_STRIDE)
    ap.add_argument("--time-stride", type=int, default=1)
    args = ap.parse_args()

    hdr = (f"{'eq':10s}{'noise':10s}{'N':>5s}{'Nt':>5s}{'m':>5s}{'p':>4s}"
           f"{'m/N':>6s}{'m_t':>5s}{'p_t':>4s}{'k_x':>6s}{'k_t':>5s}"
           f"{'cap':>6s}{'wins':>6s}{'rows':>10s}")
    print(f"{args.data_dir}   tau_hat={args.tau_hat}  tau={args.tau}")
    print(hdr)
    print("-" * len(hdr))
    for eq in args.equations:
        for label, noise in NOISE_LABELS:
            path = os.path.join(args.data_dir, f"{eq}_{label}.npz")
            if not os.path.exists(path):
                continue
            snaps = np.load(path)["snapshots"]
            _, skip_frac, _ = S._get_meta(eq)
            st = snaps[max(0, int(skip_frac * len(snaps))):]
            r = wf.mb_test_function_params(
                st, spatial_axes=range(1, st.ndim), time_axis=0,
                alpha_x=S.ALPHA_X, alpha_t=S.ALPHA_T,
                tau=args.tau, tau_hat=args.tau_hat,
            )
            Nt, N = st.shape[0], st.shape[1]
            cap = []
            if r["m"] == (min(st.shape[1:]) - 1) // 2:
                cap.append("x")
            if r["m_t"] == (Nt - 1) // 2:
                cap.append("t")
            # periodic boundary: every grid point is a centre, so the spatial
            # count is the strided grid; time keeps the valid windows only.
            n_space = int(np.prod([len(range(0, n, args.space_stride))
                                   for n in st.shape[1:]]))
            n_win = len(range(0, max(Nt - 2 * r["m_t"], 0), args.time_stride))
            print(f"{eq:10s}{label:10s}{N:5d}{Nt:5d}{r['m']:5d}{r['p']:4d}"
                  f"{r['m'] / N:6.2f}{r['m_t']:5d}{r['p_t']:4d}"
                  f"{r['k_x'][0]:6d}{r['k_t']:5d}"
                  f"{('+'.join(cap) if cap else '-'):>6s}"
                  f"{n_win:6d}{n_space * n_win:10d}")
        print()


if __name__ == "__main__":
    main()
