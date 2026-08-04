"""Minimal Modal smoke test: parallel square-root computation."""
import modal

app = modal.App("pdedisc-smoke-test")
image = modal.Image.debian_slim().pip_install("numpy")


@app.function(image=image)
def compute_sqrt(i: int) -> tuple[int, float]:
    import numpy as np
    return i, float(np.sqrt(i))


@app.local_entrypoint()
def main():
    import time

    inputs = list(range(100))
    t0 = time.time()
    results = list(compute_sqrt.map(inputs))
    elapsed = time.time() - t0

    print(f"Ran {len(inputs)} tasks in {elapsed:.1f}s")
    print(f"Sample results: {results[:5]} ... {results[-3:]}")
    print("Modal smoke test: OK")
