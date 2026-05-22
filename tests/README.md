Tests for derivative backends
=============================

This directory contains pytest suites for the two differentiation backends used by the project:

- `test_spectral_derivatives.py` — tests for the FFT-based `spectral_derivatives.py`.
- `test_polynomial_derivatives.py` — tests for the local least-squares polynomial differentiation backend `polynomial_derivatives.py`.

Running
-------
From the repository root run the following commands. These paths are repository-relative and do not depend on your current working directory.

```bash
# install dependencies for running the tests (recommended)
pip install -r requirements.txt

# run the entire tests directory
pytest -q tests

# run a single test file
pytest -q tests/test_spectral_derivatives.py
pytest -q tests/test_polynomial_derivatives.py
```

Notes
-----
- The spectral tests exercise analytic trigonometric modes and use tight tolerances (near machine epsilon) because the methods are exact for band-limited inputs.
- The polynomial backend uses local least-squares fits and is approximate; the tests therefore use relaxed tolerances and lower spatial frequencies so the local stencils approximate the analytic derivatives well.
- The tests intentionally check a few edge-cases: rectangular domains, complex-valued input handling (spectral), constant-field behaviour (both backends), and expected shape errors when a 1-D array is passed to a 2-D routine.
