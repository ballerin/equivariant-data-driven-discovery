"""
Modal-parallelised sweep: fan out by library, merge results.

Each Modal worker runs sweep.py for ONE library across the specified
equations and noise levels, writes a partial CSV, and returns the rows.
The local entrypoint merges all partial CSVs into a single result file.

Usage
-----
    # Test with Allen-Cahn, clean data, all libraries
    python scripts/sweep_modal.py --equations ac --noise 0.0

    # Full sweep, STRidge only, specific libraries
    python scripts/sweep_modal.py \\
        --equations ac sh \\
        --noise 0.0 0.02 \\
        --libraries multigraph pdefind \\
        --methods pointwise \\
        --algorithms stridge

    # Pass through any sweep.py argument
    python scripts/sweep_modal.py --equations ac --backends spectral --name test
"""

import csv
import datetime
import io
import os
import sys
import time

# ── path setup (for local imports before Modal takes over) ──────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
DATA_DIR = os.path.join(CODE_DIR, "data")
RESULTS_DIR = os.path.join(CODE_DIR, "results")

# These are the libraries sweep.py accepts
ALL_LIBRARIES = ["multigraph", "mod_multigraph", "pdefind", "polyjet", "wsindy"]

# ── Modal app definition ─────────────────────────────────────────────────

import modal

app = modal.App("pdedisc-sweep")

# Container image: Debian + numpy + the equivariant_pde_find package
# + scripts + data (mounted at container start — add_local_* MUST be last)
image = (
    modal.Image.debian_slim()
    .pip_install("numpy")
    .env({"PYTHONPATH": "/src"})
    .add_local_dir(SRC_DIR, remote_path="/src")
    .add_local_dir(os.path.join(CODE_DIR, "scripts"), remote_path="/scripts")
    .add_local_dir(DATA_DIR, remote_path="/data")
)


@app.function(
    image=image,
    timeout=7200,  # 120 min per shard
    cpu=1.0,
    memory=2048,
)
def sweep_one_shard(
    equation: str,
    library: str,
    noise: float,
    backend: str,
    extra_args: list[str],
) -> tuple[str, str, bool]:
    """Run sweep.py for one (equation, library, noise, backend) shard.

    Returns (equation, csv_content, timed_out).  timed_out=True means
    the 3500s subprocess timeout was hit — the CSV is header-only.
    """
    import subprocess

    out_path = f"/tmp/sweep_{equation}_{library}_n{noise}_{backend}.csv"

    cmd = [
        sys.executable, "/scripts/sweep.py",
        "--equations", equation,
        "--libraries", library,
        "--noise", str(noise),
        "--backends", backend,
        "--out", out_path,
        "--data-dir", "/data",
    ]
    cmd.extend(extra_args)

    label = f"[{equation}/{library}/n={noise}/{backend}]"
    print(f"{label} Running sweep.py", flush=True)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7100)
    except subprocess.TimeoutExpired:
        print(f"{label} ⚠ TIMEOUT after 7100s — skipping", flush=True)
        # Return header-only CSV so merge doesn't break
        return equation, (
            "equation,noise,method,algorithm,library,backend,weak_m,lam,"
            "d_tol,mstls_alpha,n_active,fp,fn,exact_recovery,E_coeff,"
            "max_coeff_err,E_data,cond_theta,n_samples,post_hoc_noise,"
            "spectral_cutoff,profile_type,profile_direction,time_ms\n"
        ), True

    if result.stdout:
        for line in result.stdout.strip().split("\n")[-5:]:
            print(f"  {label} {line}", flush=True)

    if result.returncode != 0:
        print(f"{label} STDERR:\n{result.stderr}", flush=True)
        raise RuntimeError(
            f"sweep.py failed for {equation}/{library}/n={noise}/{backend} "
            f"(rc={result.returncode})"
        )

    if not os.path.exists(out_path):
        raise RuntimeError(f"Output file not found: {out_path}")

    with open(out_path, "r") as f:
        csv_content = f.read()

    n_lines = csv_content.count("\n")
    print(f"{label} Done: {n_lines} lines ({n_lines - 1} data rows)",
          flush=True)
    return equation, csv_content, False


# ── local entrypoint ─────────────────────────────────────────────────────


@app.local_entrypoint()
def main(
    equations: str = "ac",
    noise: str = "0.0",
    libraries: str = "multigraph,pdefind,wsindy",
    algorithms: str = "stridge",
    methods: str = "pointwise",
    backends: str = "spectral,polynomial",
    name: str = "",
):
    """Modal-parallelised sweep — one CSV per equation.

    Parameters are passed as CLI args::

        modal run scripts/sweep_modal.py \\
            --equations ac,sh --noise 0.0,0.02 \\
            --libraries multigraph,pdefind \\
            --algorithms stridge --methods pointwise --name myrun

    Each equation gets its own CSV: results/sweep_myrun_ac_<ts>.csv, etc.
    """
    eq_list = [e.strip() for e in equations.split(",") if e.strip()]
    noise_list = [float(n.strip()) for n in noise.split(",") if n.strip()]
    lib_list = [l.strip() for l in libraries.split(",") if l.strip()]
    algo_list = [a.strip() for a in algorithms.split(",") if a.strip()]
    method_list = [m.strip() for m in methods.split(",") if m.strip()]
    backend_list = [b.strip() for b in backends.split(",") if b.strip()]

    extra_args = []
    if method_list:
        extra_args.extend(["--methods"] + method_list)
    if algo_list:
        extra_args.extend(["--algorithms"] + algo_list)

    # ── Fan-out grid: one shard per (equation, library, noise, backend) ─
    import itertools

    shards = list(itertools.product(eq_list, lib_list, noise_list, backend_list))
    eqs, libs, noises, bks = zip(*shards)

    print(f"{'='*70}")
    print(f"Modal sweep: {len(shards)} shards "
          f"({len(eq_list)} eqs × {len(lib_list)} libs × "
          f"{len(noise_list)} noise × {len(backend_list)} backends)")
    print(f"Equations: {eq_list}")
    print(f"Libraries: {lib_list}")
    print(f"Noise:     {noise_list}")
    print(f"Backends:  {backend_list}")
    print(f"Methods:   {method_list}")
    print(f"Algos:     {algo_list}")
    print(f"{'='*70}")

    t0 = time.time()

    results = list(sweep_one_shard.map(
        list(eqs), list(libs), list(noises), list(bks),
        kwargs=dict(extra_args=extra_args),
    ))
    # results is list of (equation, csv_content, timed_out) tuples

    elapsed = time.time() - t0
    print(f"\nAll {len(shards)} workers finished in {elapsed:.1f}s "
          f"({elapsed/60:.1f} min)")

    # ── Check for timeouts ──────────────────────────────────────────
    timed_out = [(eq, csv) for eq, csv, to in results if to]
    if timed_out:
        print(f"\n⚠️  WARNING: {len(timed_out)} shard(s) TIMED OUT:")
        for eq, _ in timed_out:
            print(f"     {eq}")
        print("   Missing data in the CSVs below. Re-run with higher timeout.\n")

    # ── Group by equation and merge ────────────────────────────────────
    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    from collections import defaultdict
    by_eq = defaultdict(list)
    for eq, csv_str, _to in results:
        by_eq[eq].append(csv_str)

    for eq in sorted(by_eq):
        csv_list = by_eq[eq]
        if name:
            stem = f"sweep_{name}_{eq}_{timestamp}"
        else:
            stem = f"sweep_{eq}_{timestamp}"
        out_path = os.path.join(RESULTS_DIR, f"{stem}.csv")

        total_rows = 0
        with open(out_path, "w", newline="") as fout:
            for idx, csv_str in enumerate(csv_list):
                reader = csv.reader(io.StringIO(csv_str))
                header = next(reader)
                if idx == 0:
                    writer = csv.writer(fout)
                    writer.writerow(header)
                for row in reader:
                    writer.writerow(row)
                    total_rows += 1

        print(f"  {eq}: {total_rows} rows → {os.path.basename(out_path)}")

    # Write log
    log_stem = f"sweep_{name}_{timestamp}" if name else f"sweep_{timestamp}"
    log_path = os.path.join(RESULTS_DIR, f"{log_stem}.log")
    with open(log_path, "w") as log:
        log.write(f"Modal parallel sweep\n")
        log.write(f"Timestamp: {timestamp}\n")
        log.write(f"Equations: {eq_list}\n")
        log.write(f"Libraries: {lib_list}\n")
        log.write(f"Noise levels: {noise_list}\n")
        log.write(f"Backends: {backend_list}\n")
        log.write(f"Methods: {method_list}\n")
        log.write(f"Algorithms: {algo_list}\n")
        log.write(f"Total shards: {len(shards)}\n")
        log.write(f"Wall time: {elapsed:.1f}s ({elapsed/60:.1f} min)\n")
        total_all = sum(by_eq[eq].__len__() for eq in by_eq)
        for eq in sorted(by_eq):
            log.write(f"  {eq}: {sum(1 for _ in io.StringIO(''.join(by_eq[eq][0])).readlines()) if False else 'see CSV'} rows\n")
    print(f"Log → {log_path}")
    print(f"Wall time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
