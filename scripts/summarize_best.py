"""
Summarize the best recovery from a unified sweep CSV.

Reads the sweep CSV, groups by (equation, noise, method, library, backend),
selects the best row per group using the lexicographic criterion, and
writes a concise summary CSV + prints a formatted table.

Usage
-----
    cd paper-experimental/code
    python scripts/summarize_best.py results/unified_sweep_20260607_120000.csv
    python scripts/summarize_best.py results/unified_sweep_20260607_120000.csv \\
        --group-by equation noise method library
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, CODE_DIR)

from evaluation import best_of_sweep  # noqa: E402
import math                                      # noqa: E402


def _is_nan(val):
    """Check if a value is NaN (float or string 'nan')."""
    if val is None or val == "":
        return False
    try:
        return math.isnan(float(val))
    except (ValueError, TypeError):
        return False


def summarize(csv_path, group_keys):
    """Read sweep CSV, group, select best per group.

    Parameters
    ----------
    csv_path : str
    group_keys : list of str
        Columns to group by before selecting best.

    Returns
    -------
    best_rows : list of dict
    all_rows : list of dict (all rows read)
    """
    rows = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for key in ["noise", "fp", "fn", "exact_recovery", "E",
                        "max_coeff_err", "cond_theta", "n_active",
                        "n_samples", "lam", "d_tol"]:
                if key in row and row[key] != "" and row[key] is not None:
                    try:
                        row[key] = float(row[key])
                    except (ValueError, TypeError):
                        pass
            # exact_recovery as int
            if "exact_recovery" in row:
                try:
                    row["exact_recovery"] = int(float(row["exact_recovery"]))
                except (ValueError, TypeError):
                    row["exact_recovery"] = 0
            # fp, fn, n_active as int
            for key in ["fp", "fn", "n_active"]:
                if key in row:
                    try:
                        row[key] = int(float(row[key]))
                    except (ValueError, TypeError):
                        pass
            rows.append(row)

    # Group
    groups = defaultdict(list)
    for row in rows:
        group_tuple = tuple(row.get(k, "") for k in group_keys)
        groups[group_tuple].append(row)

    # Best per group
    best_rows = []
    for group_tuple, group_rows in groups.items():
        best = best_of_sweep(group_rows, criterion="lexicographic")
        best_rows.append(best)

    # Sort for readability
    sort_keys = group_keys + ["E"]
    def sort_key(r):
        return tuple(r.get(k, "") for k in sort_keys)
    best_rows.sort(key=sort_key)

    return best_rows, rows


def print_table(best_rows, group_keys):
    """Print a human-readable summary table."""
    # Columns to display
    display_cols = group_keys + [
        "method", "algorithm", "wsindy_m",
        "d_tol", "fp", "fn", "exact_recovery", "E", "cond_theta",
    ]
    # Filter to columns that exist
    available = [c for c in display_cols if c in best_rows[0]] if best_rows else []

    # Header
    header = "  ".join(f"{c:>12}" if c not in ("exact_recovery", "fp", "fn")
                       else f"{c:>5}" for c in available)
    sep = "  ".join("-" * max(len(c), 12) for c in available)
    print(header)
    print(sep)

    for row in best_rows:
        parts = []
        for c in available:
            val = row.get(c, "")
            if c == "E":
                if _is_nan(val):
                    parts.append(f"{'nan':>12}")
                else:
                    parts.append(f"{float(val):12.4e}")
            elif c == "cond_theta":
                if _is_nan(val):
                    parts.append(f"{'nan':>12}")
                else:
                    parts.append(f"{float(val):12.2e}")
            elif c in ("exact_recovery", "fp", "fn"):
                if _is_nan(val):
                    parts.append(f"{'nan':>5}")
                elif val == "" or val is None:
                    parts.append(f"{'':>5}")
                else:
                    parts.append(f"{int(float(val)):>5}")
            elif c == "d_tol":
                if _is_nan(val):
                    parts.append(f"{'nan':>12}")
                else:
                    parts.append(f"{float(val):12.4f}")
            else:
                parts.append(f"{str(val):>12}")
        print("  ".join(parts))

    print(f"\n{len(best_rows)} groups shown.")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize best recoveries from a unified sweep CSV."
    )
    parser.add_argument(
        "csv_path",
        help="Path to the unified sweep CSV.",
    )
    parser.add_argument(
        "--group-by", nargs="+",
        default=["equation", "noise", "library"],
        help="Columns to group by (default: equation noise library).",
    )
    parser.add_argument(
        "--output", default=None,
        help="Path to save summary CSV (default: derived from input).",
    )
    parser.add_argument(
        "--table-only", action="store_true",
        help="Print table only, skip writing CSV.",
    )
    args = parser.parse_args()

    group_keys = args.group_by
    best_rows, all_rows = summarize(args.csv_path, group_keys)

    if not best_rows:
        print("No rows found.")
        return

    # Statistics
    n_total = len(all_rows)
    n_groups = len(best_rows)
    n_exact = sum(1 for r in best_rows if r.get("exact_recovery"))
    print(f"Total rows: {n_total}  Groups: {n_groups}  "
          f"Exact recovery: {n_exact}/{n_groups}")
    print()

    print_table(best_rows, group_keys)

    # Save CSV
    if not args.table_only:
        if args.output:
            out_path = args.output
        else:
            base = os.path.splitext(args.csv_path)[0]
            out_path = f"{base}_best.csv"

        if best_rows:
            fieldnames = list(best_rows[0].keys())
            with open(out_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(best_rows)
            print(f"\nSummary → {out_path}")


if __name__ == "__main__":
    main()
