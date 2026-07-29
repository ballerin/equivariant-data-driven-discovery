Tests
=====

This directory contains the pytest suite for `equivariant_pde_find`. Coverage
spans differentiation backends, library builders (multigraph, modified
multigraph, PDE-FIND, WSINDy, polynomial jet), regression algorithms, weak-form
infrastructure, PDE simulators, ground-truth equation metadata, and
end-to-end recovery pipelines (including the n-dimensional extension).

| File | Covers |
|---|---|
| `test_3d_integration.py` | End-to-end 3-D Allen-Cahn PDE discovery |
| `test_conversions_nd.py` | `library/conversions.py` — jet decompositions at arbitrary ndim, frozen 2-D regression fixture |
| `test_derivatives_backend.py` | `derivatives/__init__.py` — backend switching, dispatch, canonical `axis_lengths` |
| `test_equations.py` | `equations.py` — ground-truth coefficients and equation registry |
| `test_equivariance.py` | O(2)-equivariance of multigraph operators and the equivariant library (2-D) |
| `test_equivariance_nd.py` | O(n)-invariance of the jet decompositions (grid-free, ndim 2–6) and hyperoctahedral equivariance of the numeric operators (ndim 3–4) |
| `test_evaluation.py` | `evaluation.py` — recovery evaluation against ground truth |
| `test_jet_algebra_cas.py` | Independent **sympy** cross-check of all 40 operator jet decompositions at ndim 2–4 |
| `test_ks_integrator.py` | Kuramoto-Sivashinsky ETDRK4/IMEX integrators — coefficients, convergence order |
| `test_library.py` | `library/pdefind.py` — coordinate-basis library builders |
| `test_library_mapping.py` | `library/conversions.py` — Leibniz expansion, WSINDy-to-jet mapping |
| `test_library_nd.py` | All 5 library builders at ndim=4; numeric-column ↔ decomposition closure; published PDE-FIND/WSINDy conventions |
| `test_markdown_tables.py` | `markdown_tables.py` — table builder and shared helpers |
| `test_mod_multigraph.py` | Modified multigraph (IBP-friendly) operator basis |
| `test_mstls.py` | MSTLS — Modified Sequential Thresholding Least Squares |
| `test_multigraph_operators.py` | `multigraph_operators.py` — operator generation and ordering |
| `test_multigraph_operators_nd.py` | All 20 `max_p<=2` operators against closed-form single-mode values at ndim 4–5; `max_p=3` guard |
| `test_ndim.py` | N-dimensional PDE discovery extensions |
| `test_polynomial_derivatives.py` | Local least-squares polynomial derivative backend, incl. N-D convergence order |
| `test_regression.py` | `regression.py` — STRidge, TrainSTRidge |
| `test_rfe.py` | RFE — Recursive Feature Elimination |
| `test_simulations.py` | PDE simulators — smoke tests and physical correctness |
| `test_spectral_derivatives.py` | FFT-based spectral derivative backend, incl. Nyquist convention and N-D exactness |
| `test_utils.py` | `utils.py` — `print_pde`, `rotation_90`, anisotropic noise |
| `test_weak_form.py` | `library/weak_form.py` — weak-form PDE discovery |
| `test_weak_form_nd.py` | N-D weak form — exact IBP kernel expansions, analytic single-mode checks, `compute_m_scale` |
| `test_wsindy_jet_libraries.py` | WSINDy library and polynomial jet library |

How the N-D generalization is verified
--------------------------------------

The operator basis is defined once symbolically (`jet_algebra.py`) and that one
definition drives both the numeric library column and the analytic jet
decomposition. That rules out silent divergence between them, but it also means
most cross-checks compare the code against itself. The tests above therefore
supply **independent** ground truth at ndim > 2 from four directions:

1. **An outside CAS** — `test_jet_algebra_cas.py` re-derives every operator with
   sympy, importing nothing from `jet_algebra`.
2. **The defining mathematical property** — `test_equivariance_nd.py` checks
   O(n)-invariance directly under random `SO(n)` rotations, which rejects
   formulas that lattice symmetry alone would accept.
3. **Closed-form analytic values** — `test_multigraph_operators_nd.py` compares
   against exact single-Fourier-mode derivatives with *distinct* per-axis
   wavenumbers, so an axis-ordering error cannot hide behind the isotropic
   domain used beyond 3-D.
4. **Independently written code paths agreeing** —
   `test_library_nd.py::TestDecompositionMatchesNumericLibrary` forces the
   numeric library builder, the decomposition table and the jet library builder
   to produce the same column at ndim = 2, 3, 4.

Several of these classes contain an explicit "is this check discriminating?"
test that feeds the comparison a deliberately wrong input and requires it to
fail — a cross-check that cannot fail proves nothing.

Running
-------
From the repository root:

```bash
# install the package with test dependencies
pip install -e ".[test]"

# run the entire test suite
pytest -q tests

# run a single test file
pytest -q tests/test_spectral_derivatives.py

# skip slow/integration tests
pytest -q tests -m "not slow and not integration"
```

Notes
-----
- The spectral derivative tests exercise analytic trigonometric modes and use tight tolerances (near machine epsilon) because the methods are exact for band-limited inputs.
- The polynomial backend uses local least-squares fits and is approximate; the tests therefore use relaxed tolerances and lower spatial frequencies so the local stencils approximate the analytic derivatives well.
- A few tests in `test_mstls.py` are marked `integration`/`slow` — they run full simulate → library → regression → recovery pipelines and take noticeably longer than the rest of the suite.
