"""
Modal-parallelised N-D sweep: fan out by (equation, ndim, noise, library).

Each Modal worker runs sweep_nd.py for ONE (equation, ndim, noise, library)
combination with the appropriate memory allocation.  The local entrypoint
merges all partial CSVs into a single result file.

Resources are allocated per (ndim, library) based on measured peak RSS
from local 2D-4D runs and conservative extrapolation for 5D.

Usage
-----
    # Local test (runs on Modal, but fan-out is 1 shard)
    python scripts/sweep_nd_modal.py --equations ac --ndim 2 --noise 0.0 --libraries multigraph

    # Full N-D sweep on Modal
    modal run scripts/sweep_nd_modal.py --name ndim_full

    # 5D only (heavy)
    modal run scripts/sweep_nd_modal.py --ndim 5 --name ndim_5d
"""

import csv
import datetime
import io
import os
import sys
import time

# ── path setup ──────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
SRC_DIR = os.path.join(CODE_DIR, "src")
DATA_DIR = os.path.join(CODE_DIR, "data")
RESULTS_DIR = os.path.join(CODE_DIR, "results")

ALL_LIBRARIES = ["multigraph", "mod_multigraph", "pdefind", "wsindy"]

# ── Resource allocation table ───────────────────────────────────────────
# (ndim, library) → (memory_mb, timeout_s)
# Based on measured peak RSS from local 2D-4D runs (2026-08-04).
# 5D values are conservative extrapolations (2× 4D peak for multigraph,
# 1.5× for coordinate libraries since column count grows sub-linearly
# with grid size at fixed max_order=4).
RESOURCES = {
    # 2D: everything fits in 256 MB
    (2, "multigraph"):     (256,   600),
    (2, "mod_multigraph"): (256,   600),
    (2, "pdefind"):        (256,   600),
    (2, "wsindy"):         (256,   600),
    # 3D: peak ~283 MB → 512 MB with headroom
    (3, "multigraph"):     (512,  1200),
    (3, "mod_multigraph"): (512,  1200),
    (3, "pdefind"):        (512,  1200),
    (3, "wsindy"):         (512,  1200),
    # 4D: multigraph ~1.2 GB → 2048; pdefind ~2.3 GB → 4096; wsindy ~3.2 GB → 4096
    (4, "multigraph"):     (2048,  5400),
    (4, "mod_multigraph"): (2048,  5400),
    (4, "pdefind"):        (4096,  5400),
    (4, "wsindy"):         (4096,  5400),
    # 5D: conservative extrapolation
    (5, "multigraph"):     (4096, 14400),
    (5, "mod_multigraph"): (4096, 14400),
    (5, "pdefind"):        (8192, 14400),
    (5, "wsindy"):         (10240, 14400),
}

# ── Modal app ───────────────────────────────────────────────────────────

import modal

app = modal.App("pdedisc-sweep-nd")

# Persistent volume for sweep results (survives local disconnects)
volume = modal.Volume.from_name("pdedisc-sweep-results", create_if_missing=True)

image = (
    modal.Image.debian_slim()
    .pip_install("numpy")
    .env({"PYTHONPATH": "/src", "OMP_NUM_THREADS": "1"})
    .add_local_dir(SRC_DIR, remote_path="/src")
    .add_local_dir(os.path.join(CODE_DIR, "scripts"), remote_path="/scripts")
    .add_local_dir(DATA_DIR, remote_path="/data")
)


def _run_shard(equation, ndim, noise, library, extra_args, label, timeout_s):
    """Shared worker logic — runs sweep_nd.py, returns (eq, csv, timed_out)."""
    import subprocess

    out_path = f"/tmp/sweep_nd_{equation}_nd{ndim}_n{noise}_{library}.csv"
    cmd = [
        sys.executable, "/scripts/sweep_nd.py",
        "--equations", equation,
        "--ndim", str(ndim),
        "--noise", str(noise),
        "--libraries", library,
        "--out", out_path,
        "--data-dir", "/data",
    ]
    cmd.extend(extra_args)

    tag = f"[{equation}/nd{ndim}/n={noise}/{library}/{label}]"
    print(f"{tag} Running sweep_nd.py", flush=True)
    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=max(10, timeout_s - 60),
        )
    except subprocess.TimeoutExpired:
        dt = time.time() - t0
        print(f"{tag} TIMEOUT after {dt:.0f}s", flush=True)
        return equation, (
            "equation,ndim,noise,library,algorithm,use_ibp,m,d_tol,lam,"
            "n_active,fp,fn,exact_recovery,E_coeff,E_data,cond_theta,"
            "n_cols,n_rows,n_snapshots,grid_points,"
            "build_time_ms,fit_time_ms,peak_rss_mb\n"
        ), True

    dt = time.time() - t0
    if result.stdout:
        for line in result.stdout.strip().split("\n")[-5:]:
            print(f"  {tag} {line}", flush=True)
    if result.returncode != 0:
        print(f"{tag} STDERR:\n{result.stderr}", flush=True)
        raise RuntimeError(
            f"sweep_nd.py failed for {equation}/nd{ndim}/n={noise}/{library}"
        )
    if not os.path.exists(out_path):
        raise RuntimeError(f"Output not found: {out_path}")

    with open(out_path, "r") as f:
        csv_content = f.read()

    # Persist CSV to Modal Volume so results survive local disconnects
    vol_path = f"/results/sweep_nd_{equation}_nd{ndim}_n{noise}_{library}.csv"
    try:
        with open(vol_path, "w") as vf:
            vf.write(csv_content)
    except Exception as e:
        print(f"  {tag} WARNING: failed to write to Volume: {e}", flush=True)

    n_lines = csv_content.count("\n")
    print(f"{tag} Done: {n_lines} lines in {dt:.0f}s", flush=True)
    return equation, csv_content, False


# ── Modal functions at global scope (one per memory tier) ──────────────

@app.function(image=image, timeout=600,  cpu=1.0, memory=256,  volumes={"/results": volume})
def w256(equation, ndim, noise, library, extra_args):
    return _run_shard(equation, ndim, noise, library, extra_args, "m256", 600)

@app.function(image=image, timeout=1200, cpu=1.0, memory=512,  volumes={"/results": volume})
def w512(equation, ndim, noise, library, extra_args):
    return _run_shard(equation, ndim, noise, library, extra_args, "m512", 1200)

@app.function(image=image, timeout=5400, cpu=1.0, memory=2048, volumes={"/results": volume})
def w2048(equation, ndim, noise, library, extra_args):
    return _run_shard(equation, ndim, noise, library, extra_args, "m2048", 5400)

@app.function(image=image, timeout=5400, cpu=1.0, memory=4096, volumes={"/results": volume})
def w4096(equation, ndim, noise, library, extra_args):
    return _run_shard(equation, ndim, noise, library, extra_args, "m4096", 5400)


def _get_worker(ndim, library):
    mem, _ = RESOURCES.get((ndim, library), (4096, 5400))
    if mem <= 256:    return w256
    elif mem <= 512:  return w512
    elif mem <= 2048: return w2048
    else:             return w4096


# ── local entrypoint ────────────────────────────────────────────────────

@app.local_entrypoint()
def main(
    equations: str = "ac,sh,pme,ch",
    ndim: str = "2,3,4,5",
    noise: str = "0.0,0.02,0.05",
    libraries: str = "multigraph,mod_multigraph,pdefind,wsindy",
    algorithms: str = "stride,rfe",
    ibp: str = "both",
    name: str = "",
):
    """Modal-parallelised N-D sweep.

    Parameters are passed as comma-separated strings (Modal convention)::

        modal run scripts/sweep_nd_modal.py --equations ac --ndim 2,3,4,5
        modal run scripts/sweep_nd_modal.py --ndim 5 --name ndim_5d_only
    """
    eq_list = [e.strip() for e in equations.split(",") if e.strip()]
    ndim_list = [int(n.strip()) for n in ndim.split(",") if n.strip()]
    noise_list = [float(n.strip()) for n in noise.split(",") if n.strip()]
    lib_list = [l.strip() for l in libraries.split(",") if l.strip()]
    algo_list = [a.strip() for a in algorithms.split(",") if a.strip()]

    extra_args = []
    if algo_list:
        extra_args.extend(["--algorithms"] + algo_list)
    if ibp:
        extra_args.extend(["--ibp", ibp])

    # ── Fan-out grid ──────────────────────────────────────────────────
    import itertools

    shards = list(itertools.product(eq_list, ndim_list, noise_list, lib_list))
    # PDE-FIND cannot represent PME (products of derivatives) — skip.
    shards = [(e, n, no, l) for e, n, no, l in shards
              if not (e == "pme" and l == "pdefind")]
    eqs, ndims, noises, libs = zip(*shards)

    # ── Resource summary ──────────────────────────────────────────────
    total_mem = sum(
        RESOURCES.get((nd, lb), (4096, 14400))[0] / 1024
        for nd, lb in zip(ndims, libs)
    )
    max_timeout = max(
        RESOURCES.get((nd, lb), (4096, 14400))[1]
        for nd, lb in zip(ndims, libs)
    )

    print(f"{'='*70}")
    print(f"Modal N-D sweep: {len(shards)} shards")
    print(f"  Equations: {eq_list}")
    print(f"  Dimensions: {ndim_list}")
    print(f"  Noise: {noise_list}")
    print(f"  Libraries: {lib_list}")
    print(f"  Algorithms: {algo_list}")
    print(f"  IBP: {ibp}")
    print(f"  Total memory if all parallel: {total_mem:.1f} GB")
    print(f"  Max timeout: {max_timeout}s ({max_timeout/3600:.1f}h)")
    print(f"{'='*70}")

    # ── Launch shards (grouped by worker tier) ────────────────────────
    t0 = time.time()

    # Group shards by worker function so we can call .map() per tier
    from collections import defaultdict
    by_worker = defaultdict(list)
    for i, (eq, nd, noise_val, lib) in enumerate(shards):
        worker = _get_worker(nd, lib)
        by_worker[worker].append(i)

    all_results = [None] * len(shards)

    for worker, indices in by_worker.items():
        tier_eqs = [eqs[i] for i in indices]
        tier_ndims = [ndims[i] for i in indices]
        tier_noises = [noises[i] for i in indices]
        tier_libs = [libs[i] for i in indices]

        print(f"\n  [{len(indices)} shards] Launching...", flush=True)
        tier_results = list(worker.map(
            tier_eqs, tier_ndims, tier_noises, tier_libs,
            kwargs=dict(extra_args=extra_args),
        ))

        for j, idx in enumerate(indices):
            all_results[idx] = tier_results[j]

    elapsed = time.time() - t0
    print(f"\nAll {len(shards)} workers finished in {elapsed:.1f}s "
          f"({elapsed/60:.1f} min)")

    # ── Check for timeouts ────────────────────────────────────────────
    timed_out = [
        (eqs[i], ndims[i], libs[i])
        for i, (_, _, to) in enumerate(all_results) if to
    ]
    if timed_out:
        print(f"\n  {len(timed_out)} shard(s) TIMED OUT:")
        for eq, nd, lib in timed_out:
            print(f"     {eq}/nd{nd}/{lib}")

    # ── Merge all CSVs ────────────────────────────────────────────────
    os.makedirs(RESULTS_DIR, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = f"sweep_nd_{name}_{timestamp}" if name else f"sweep_nd_{timestamp}"
    out_path = os.path.join(RESULTS_DIR, f"{stem}.csv")

    total_rows = 0
    with open(out_path, "w", newline="") as fout:
        for idx, (eq, csv_str, _to) in enumerate(all_results):
            reader = csv.reader(io.StringIO(csv_str))
            try:
                header = next(reader)
            except StopIteration:
                continue
            if idx == 0 or total_rows == 0:
                writer = csv.writer(fout)
                writer.writerow(header)
            for row in reader:
                if row:  # skip empty lines
                    writer.writerow(row)
                    total_rows += 1

    print(f"\n  Merged: {total_rows} rows → {os.path.basename(out_path)}")

    # ── Write log ─────────────────────────────────────────────────────
    log_path = out_path.replace(".csv", ".log")
    with open(log_path, "w") as log:
        log.write(f"Modal N-D sweep\n")
        log.write(f"Timestamp: {timestamp}\n")
        log.write(f"Equations: {eq_list}\n")
        log.write(f"Dimensions: {ndim_list}\n")
        log.write(f"Noise levels: {noise_list}\n")
        log.write(f"Libraries: {lib_list}\n")
        log.write(f"Algorithms: {algo_list}\n")
        log.write(f"IBP: {ibp}\n")
        log.write(f"Total shards: {len(shards)}\n")
        log.write(f"Wall time: {elapsed:.1f}s ({elapsed/60:.1f} min)\n")
        for (nd, lib), (mem, timeout) in sorted(RESOURCES.items()):
            log.write(f"  ndim={nd} {lib}: memory={mem}MB timeout={timeout}s\n")
    print(f"  Log → {os.path.basename(log_path)}")
    print(f"  Wall time: {elapsed:.1f}s ({elapsed/60:.1f} min)")


# ── Merge-only function (callable independently after detached run) ─────

@app.function(image=image, timeout=600, cpu=1.0, memory=1024,
              volumes={"/results": volume})
def merge_results(name: str = "", results_dir: str = "/results"):
    """Merge all per-shard CSVs from the Volume into one result file."""
    import datetime
    import glob

    os.makedirs(RESULTS_DIR_LOCAL := "/results", exist_ok=True)

    csv_files = sorted(glob.glob(f"{results_dir}/sweep_nd_*.csv"))
    if not csv_files:
        print("No CSV files found on Volume. The sweep may not have completed.")
        return

    print(f"Found {len(csv_files)} per-shard CSV files on Volume")

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = f"sweep_nd_{name}_{timestamp}" if name else f"sweep_nd_{timestamp}"
    merged_path = f"{results_dir}/{stem}.csv"

    total_rows = 0
    with open(merged_path, "w", newline="") as fout:
        writer = None
        for i, fpath in enumerate(csv_files):
            with open(fpath, "r") as fin:
                reader = csv.reader(fin)
                try:
                    header = next(reader)
                except StopIteration:
                    continue
                if writer is None:
                    writer = csv.writer(fout)
                    writer.writerow(header)
                for row in reader:
                    if row:
                        writer.writerow(row)
                        total_rows += 1

    print(f"Merged {total_rows} rows from {len(csv_files)} files → {os.path.basename(merged_path)}")

    # Write log
    log_path = merged_path.replace(".csv", ".log")
    with open(log_path, "w") as log:
        log.write(f"Modal N-D sweep merge\n")
        log.write(f"Timestamp: {timestamp}\n")
        log.write(f"Name: {name}\n")
        log.write(f"Shard files: {len(csv_files)}\n")
        log.write(f"Total rows: {total_rows}\n")
        log.write(f"Source files:\n")
        for fpath in csv_files:
            log.write(f"  {os.path.basename(fpath)}\n")
    print(f"Log → {os.path.basename(log_path)}")

    # Commit volume changes
    try:
        volume.commit()
        print("Volume committed.")
    except Exception as e:
        print(f"Volume commit note: {e}")

    return merged_path
