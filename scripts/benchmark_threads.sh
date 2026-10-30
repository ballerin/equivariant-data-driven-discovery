#!/usr/bin/env bash
# benchmark_threads.sh — Measure sweep.py wall time vs. thread count.
#
# Usage:
#   chmod +x scripts/benchmark_threads.sh
#   bash scripts/benchmark_threads.sh
#
# The Allen–Cahn clean dataset is swept with:
#   - spectral differentiation (pointwise only)
#   - multigraph library
#   - STRidge regression (4 λ × 17 d_tol = 68 combos)
#
# Runs are repeated for 1, 2, 4, 8 threads.

set -uo pipefail  # no -e: we handle errors manually

THREAD_COUNTS=(1 2 4 8)
RESULTS_DIR="results/benchmark_threads"
mkdir -p "$RESULTS_DIR"

SUMMARY_FILE="$RESULTS_DIR/summary.txt"
echo "Thread benchmark — $(date)" > "$SUMMARY_FILE"
echo "Command: sweep.py --equations ac --noise 0.0 --methods pointwise --algorithms stridge --backends spectral --libraries multigraph --no-oracle" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
printf "%-10s %12s %12s\n" "Threads" "Wall(s)" "Rows/s" >> "$SUMMARY_FILE"
printf "%-10s %12s %12s\n" "-------" "--------" "------" >> "$SUMMARY_FILE"

for N in "${THREAD_COUNTS[@]}"; do
    echo ""
    echo "============================================================"
    echo "  Running with OPENBLAS_NUM_THREADS=$N OMP_NUM_THREADS=$N"
    echo "============================================================"

    export OPENBLAS_NUM_THREADS=$N
    export OMP_NUM_THREADS=$N
    export NUMEXPR_NUM_THREADS=$N
    export MKL_NUM_THREADS=$N
    export VECLIB_MAXIMUM_THREADS=$N

    START_TIME=$(date +%s.%N)

    python scripts/sweep.py \
        --equations ac \
        --noise 0.0 \
        --methods pointwise \
        --algorithms stridge \\
        --backends spectral \
        --libraries multigraph \
        --no-oracle \
        --name "threads_${N}"

    END_TIME=$(date +%s.%N)
    ELAPSED=$(python3 -c "print(f'{${END_TIME} - ${START_TIME}:.1f}')")

    # Count how many rows were written (excl. header)
    CSV_FILE=$(ls -t results/sweep_*_threads_${N}.csv 2>/dev/null | head -1)
    if [[ -f "$CSV_FILE" ]]; then
        NROWS=$(($(wc -l < "$CSV_FILE") - 1))
        RATE=$(python3 -c "print(f'{${NROWS} / ${ELAPSED}:.1f}')")
    else
        NROWS="?"
        RATE="?"
    fi

    printf "%-10s %12s %12s\n" "$N" "$ELAPSED" "$RATE" | tee -a "$SUMMARY_FILE"

    # Move CSV to benchmark dir
    if [[ -f "$CSV_FILE" ]]; then
        mv "$CSV_FILE" "$RESULTS_DIR/"
        LOG_FILE="${CSV_FILE%.csv}.log"
        mv "$LOG_FILE" "$RESULTS_DIR/" 2>/dev/null || true
    fi

    echo "  OK: ${ELAPSED}s  rows=${NROWS}"
done

echo ""
echo "============================================================"
echo "  Summary"
echo "============================================================"
cat "$SUMMARY_FILE"
