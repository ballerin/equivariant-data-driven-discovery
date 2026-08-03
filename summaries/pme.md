# Porous Medium equation ($m = 3$)

The equation under consideration is

$$f_t = \nabla\cdot(f^2\nabla f) = f^2\Delta f + 2f\|\nabla f\|^2 = \tfrac13\Delta(f^3)$$

which can be written in terms of $O(2)$-equivariant terms

| Operator               | True coefficient |
|:-----------------------|:-----------------|
| $f^2\cdot\Delta f$     | $+1.0$           |
| $f\cdot\|\nabla f\|^2$ | $+2.0$           |

or standard terms (PolyJet basis)

| Operator                  | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f^2\cdot\partial_{xx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{yy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_x^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |
| $f\cdot f_y^2$            |            | 🔴        | 🟢      | 🟢        | 🟢          |

In the standard convention $u_t = \Delta(u^m)$ this is $m = 3$; the $1/3$ prefactor is inessential and absorbed by rescaling time.  The coefficients below belong to the $\nabla\cdot(f^2\nabla f)$ normalisation actually simulated.  PDE-FIND cannot represent $f\cdot f_x^2$ or $f\cdot f_y^2$ (products of derivatives), so its minimum fn is 2.  WSINDy represents PME exactly via the chain rule: $f_t = \frac13\partial_x^2(f^3) + \frac13\partial_y^2(f^3)$, which expands to all 4 jet terms.

### Legend for summary table

| Symbol | Meaning |
|---|---|
| 🟢 | Equivariant **clear** win (fewer total fp+fn) |
| 🟩 | Equivariant **marginal** win (same fp+fn, lower E) |
| ⚪ | Tie (same fp+fn, same E) |
| 🟥 | Non-equivariant **marginal** win (same fp+fn, lower E) |
| 🔴 | Non-equivariant **clear** win (fewer total fp+fn) |
| — | Both methods failed, no comparison |
| ✓ | Exact recovery (fp = fn = 0) |
| ✗ | Failed |


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best       | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:---------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0011     | multigraph/RFE     | 0.0019    | WSINDy/RFE     | 0/0 ✓     | 0.0010      | 0.0019     | —            |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0008     | multigraph/RFE     | 0.0051    | PolyJet/RFE    | 0/0 ✓     | 0.0008      | 0.0051     | —            |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0008     | multigraph/RFE     | 0.0006    | WSINDy/RFE     | 0/0 ✓     | 0.0007      | 0.0007     | —            |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0007     | multigraph/RFE     | 0.0007    | WSINDy/RFE     | 0/0 ✓     | 0.0007      | 0.0007     | —            |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0009     | multigraph/RFE     | 0.0010    | WSINDy/RFE     | 0/0 ✓     | 0.0008      | 0.0011     | —            |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0009     | multigraph/STRidge | 0.0019    | PolyJet/MSTLS  | 0/0 ✓     | 0.0009      | 0.0011     | —            |
| 2%    | pointwise     | spectral   | —        | 10/0 ✗            | 0.2595     | —                  | 0.8545    | —              | 12/0 ✗    | 0.2552      | 0.8601     | —            |
| 2%    | pointwise     | polynomial | —        | 13/0 ✗            | 0.1117     | —                  | 0.8506    | —              | 3/0 ✗     | 0.0388      | 0.8533     | —            |
| 2%    | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0026     | multigraph/STRidge | 0.2406    | WSINDy/RFE     | 0/0 ✓     | 0.0017      | 0.2405     | —            |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0021     | multigraph/STRidge | 0.3302    | WSINDy/RFE     | 0/0 ✓     | 0.0017      | 0.2405     | —            |
| 2%    | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0037     | multigraph/STRidge | 0.2795    | WSINDy/STRidge | 0/0 ✓     | 0.0019      | 0.2795     | —            |
| 2%    | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0012     | multigraph/MSTLS   | 0.4409    | WSINDy/RFE     | 0/0 ✓     | 0.0010      | 0.4409     | —            |
| 5%    | pointwise     | spectral   | —        | 32/0 ✗            | 0.7804     | —                  | 0.9699    | —              | 22/0 ✗    | 0.8836      | 0.9702     | —            |
| 5%    | pointwise     | polynomial | —        | 4/0 ✗             | 0.1883     | —                  | 0.9716    | —              | 7/0 ✗     | 0.2165      | 0.9722     | —            |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0052     | multigraph/STRidge | 0.6476    | WSINDy/STRidge | 0/0 ✓     | 0.0056      | 0.5225     | —            |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0020     | multigraph/STRidge | 0.5225    | WSINDy/STRidge | 0/0 ✓     | 0.0056      | 0.5225     | —            |
| 5%    | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0107     | multigraph/MSTLS   | 0.5848    | WSINDy/STRidge | 0/0 ✓     | 0.0049      | 0.5842     | —            |
| 5%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0017     | multigraph/STRidge | 0.5842    | WSINDy/STRidge | 0/0 ✓     | 0.0048      | 0.5842     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best   |
|:----------------------|:--------------------------------|:------|--------:|---------:|-----------:|
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0007  | 0.0      | —          |
| WSINDy + STRidge      | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0007  | 0.2      | +2.8%      |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.0      | +3.2%      |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.7      | +3.2%      |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +3.2%      |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.0      | +3.2%      |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +6.3%      |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.5      | +6.3%      |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.8360  | 0.2      | +118725.8% |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 2/2 ✗ | 0.3963  | 0.3      | ✗          |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 6/2 ✗ | 0.1893  | 0.6      | ✗          |
| PDE-FIND + RFE        | spectral, pointwise             | 2/4 ✗ | 0.2585  | 0.0      | ✗          |

### 2% noise

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:----------------------|:--------------------------------|:------|--------:|---------:|----------:|
| WSINDy + RFE          | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0010  | 0.0      | —         |
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0012  | 0.1      | +31.3%    |
| Multigraphs + STRidge | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0017  | 0.1      | +82.5%    |
| WSINDy + STRidge      | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0019  | 0.2      | +102.2%   |
| PolyJet + MSTLS       | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0022  | 1.4      | +132.5%   |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0028  | 0.0      | +198.9%   |
| PolyJet + STRidge     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0042  | 0.5      | +340.4%   |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=4    | 0/0 ✓ | 0.5144  | 0.5      | +53955.3% |
| PolyJet + RFE         | spectral, weak (IBP), m=8       | 0/1 ✗ | 0.2944  | 0.0      | ✗         |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 1/2 ✗ | 0.4002  | 0.2      | ✗         |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 5/2 ✗ | 0.2326  | 0.8      | ✗         |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=4    | 2/4 ✗ | 0.2586  | 0.0      | ✗         |

### 5% noise

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best |
|:----------------------|:--------------------------------|:------|--------:|---------:|---------:|
| Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0017  | 0.1      | —        |
| WSINDy + STRidge      | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0048  | 0.2      | +183.9%  |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0056  | 0.1      | +226.2%  |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0056  | 0.0      | +226.2%  |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0059  | 1.5      | +248.6%  |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0193  | 0.6      | +1030.8% |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=4  | 0/0 ✓ | 0.0441  | 0.5      | +2486.9% |
| WSINDy + RFE          | polynomial, weak (IBP), m=12    | 0/2 ✗ | 0.5773  | 0.0      | ✗        |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=4  | 2/2 ✗ | 0.5205  | 0.4      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (no IBP), m=4  | 4/2 ✗ | 0.2377  | 0.9      | ✗        |
| PolyJet + RFE         | polynomial, weak (IBP), m=8     | 0/3 ✗ | 0.6655  | 0.1      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise             | 1/4 ✗ | 0.7000  | 0.0      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 10/2  | 12       | 0.6329  | 0.2109 | 0.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/4   | 2        | 0.2585  | 0.2565 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 24/2  | 26       | 0.4869  | 0.8578 | 0.7      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 24/2  | 26       | 1.0370  | 0.9698 | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 0/0   | 2        | 0.0011  | 0.0019 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 2        | 0.0010  | 0.0019 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.001   | 12/0  | 8        | 0.2552  | 0.8601 | 0.4      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/4   | 2        | 0.3045  | 0.8647 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 22/0  | 15       | 0.8836  | 0.9702 | 0.4      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 6/0   | 10       | 0.0016  | 0.0016 | 0.8      | 🔴      |
| PolyJet     | Clean | RFE       |         | 0/0   | 4        | 0.0011  | 0.0019 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.001   | 29/0  | 33       | 0.2409  | 0.8508 | 1.8      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 73/0  | 77       | 0.9578  | 0.9662 | 2.6      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 0/0   | 2        | 0.0011  | 0.0019 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 2        | 0.0011  | 0.0019 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.001   | 10/0  | 8        | 0.2595  | 0.8545 | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.001   | 32/0  | 17       | 0.7804  | 0.9699 | 0.3      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.001   | 16/2  | 18       | 0.6129  | 0.1748 | 0.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/4   | 2        | 0.2590  | 0.2565 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 5/2   | 7        | 0.2541  | 0.8598 | 0.5      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 27/2  | 29       | 0.7600  | 0.9706 | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 0/0   | 2        | 0.0013  | 0.0057 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 2        | 0.0012  | 0.0057 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.01    | 3/0   | 3        | 0.0388  | 0.8533 | 0.3      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 0/2   | 1        | 0.5784  | 0.9065 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.001   | 7/0   | 7        | 0.2165  | 0.9722 | 0.5      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 13/0  | 17       | 0.0029  | 0.0043 | 1.0      | 🔴      |
| PolyJet     | Clean | RFE       |         | 0/0   | 4        | 0.0008  | 0.0051 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 20/0  | 24       | 0.1182  | 0.8488 | 1.0      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 37/0  | 41       | 0.2905  | 0.9686 | 1.5      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 0/0   | 2        | 0.0009  | 0.0051 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 2        | 0.0008  | 0.0051 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.01    | 13/0  | 8        | 0.1117  | 0.8506 | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 4/0   | 5        | 0.1883  | 0.9716 | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | no  | 2/2   | 4        | 0.3963  | 0.3115 | 0.3      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 6/2   | 8        | 0.1893  | 0.1749 | 0.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 2/4   | 2        | 0.2585  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.01    | 8  | yes | 5/2   | 7        | 0.2914  | 0.2587 | 0.4      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 6/2   | 8        | 0.1926  | 0.1702 | 0.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/4   | 2        | 0.2586  | 0.2485 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | no  | 3/2   | 5        | 0.5812  | 0.5944 | 0.3      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 4  | no  | 5/2   | 7        | 0.2326  | 0.5876 | 0.8      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 2/4   | 2        | 0.2586  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 8  | yes | 4/2   | 6        | 0.2875  | 0.5053 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 6/2   | 8        | 0.2323  | 0.4824 | 0.7      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2586  | 0.4936 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 8  | no  | 4/2   | 6        | 0.2576  | 0.6662 | 0.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | no  | 4/2   | 6        | 0.2617  | 0.6660 | 0.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/4   | 1        | 0.7004  | 0.9255 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 4/2   | 6        | 0.2664  | 0.6074 | 0.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | yes | 4/2   | 6        | 0.2467  | 0.7773 | 0.6      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2587  | 0.7809 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.8369  | 0.8354 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0008  | 0.0011 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0022  | 0.2406 | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.5144  | 0.7166 | 0.5      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0017  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0019  | 0.2795 | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.6468  | 0.7323 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0021  | 0.2795 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 2        | 0.4008  | 0.7141 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 8  | no  | 0/2   | 1        | 0.5774  | 0.7821 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0049  | 0.5842 | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.1728  | 0.7787 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 0/2   | 1        | 0.5774  | 0.7483 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.4      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0010  | 0.0010 | 0.4      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0009  | 0.0010 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0009  | 0.0010 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0081  | 0.3304 | 0.6      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0042  | 0.3303 | 1.5      | 🟢      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/1   | 3        | 0.2979  | 0.4159 | 0.1      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.1     | 8  | yes | 3/0   | 7        | 0.0212  | 0.4420 | 0.6      | 🔴      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0041  | 0.2795 | 1.0      | 🟢      |
| PolyJet     | 2%    | RFE       |         | 8  | yes | 0/1   | 3        | 0.2944  | 0.5141 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 10.0    | 12 | no  | 6/0   | 10       | 0.0564  | 0.5228 | 0.5      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0109  | 0.6475 | 1.5      | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | no  | 2/3   | 3        | 0.6558  | 0.7893 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 8  | yes | 6/0   | 10       | 0.0943  | 0.7676 | 0.7      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0109  | 0.5848 | 1.0      | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/3   | 1        | 0.6806  | 0.8897 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 8  | no  | 0/0   | 2        | 0.0009  | 0.0011 | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 2        | 0.0010  | 0.0010 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0009  | 0.0010 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0009  | 0.0010 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0026  | 0.2406 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0038  | 0.2406 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0038  | 0.2406 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0037  | 0.2795 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0041  | 0.2795 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0041  | 0.2795 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 2        | 0.0052  | 0.6476 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0107  | 0.5229 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0107  | 0.5229 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 2        | 0.0117  | 0.5848 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0107  | 0.5848 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 0/0   | 2        | 0.0255  | 0.7696 | 0.0      | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 4  | no  | 3/2   | 5        | 0.3673  | 0.3697 | 0.4      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 6/2   | 8        | 0.1951  | 0.1769 | 0.6      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 2/4   | 2        | 0.2590  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | yes | 4/2   | 6        | 0.2358  | 0.2198 | 0.4      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 6/2   | 8        | 0.1994  | 0.1715 | 0.5      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/4   | 2        | 0.2591  | 0.2485 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 10.0    | 12 | no  | 1/2   | 3        | 0.4002  | 0.4377 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 4  | no  | 6/2   | 8        | 0.2308  | 0.5872 | 0.9      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 2/4   | 2        | 0.2590  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 8  | yes | 5/2   | 7        | 0.2919  | 0.4979 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 6/2   | 8        | 0.2322  | 0.4824 | 0.6      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2591  | 0.4936 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 4  | no  | 2/2   | 4        | 0.5205  | 0.8937 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 4  | no  | 4/2   | 6        | 0.2377  | 0.8569 | 0.9      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/4   | 1        | 0.7005  | 0.9255 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 8  | yes | 5/2   | 7        | 0.2937  | 0.7824 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 8  | yes | 4/2   | 6        | 0.2378  | 0.7773 | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2592  | 0.7809 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 2        | 0.0009  | 0.0012 | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.8369  | 0.8353 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0022  | 0.2406 | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.5209  | 0.7119 | 0.5      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0017  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0019  | 0.2795 | 0.2      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.6483  | 0.7311 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 8  | yes | 0/0   | 2        | 0.0010  | 0.4409 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.0441  | 0.8533 | 0.5      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 8  | no  | 0/2   | 1        | 0.5773  | 0.7821 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0048  | 0.5842 | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.1780  | 0.7766 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 12 | yes | 0/2   | 1        | 0.5773  | 0.7483 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.5      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0009  | 0.0011 | 0.5      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0009  | 0.0011 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0009  | 0.0011 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0068  | 0.2406 | 0.5      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0034  | 0.3302 | 1.5      | 🟢      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/2   | 2        | 0.3503  | 0.4647 | 0.0      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0042  | 0.2795 | 0.5      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0022  | 0.4409 | 1.4      | 🟢      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/1   | 3        | 0.2983  | 0.3767 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0193  | 0.6471 | 0.6      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0059  | 0.6469 | 1.5      | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | no  | 0/3   | 1        | 0.6700  | 0.8245 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.0332  | 0.7671 | 0.6      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0103  | 0.7668 | 1.4      | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/3   | 1        | 0.6655  | 0.8780 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 0.01    | 8  | yes | 0/0   | 2        | 0.0009  | 0.0019 | 0.2      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 2        | 0.0021  | 0.3302 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0028  | 0.2405 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0028  | 0.2405 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | yes | 0/0   | 2        | 0.0017  | 0.4409 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.0012  | 0.4409 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 8  | yes | 3/4   | 1        | 0.6830  | 0.7657 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 2        | 0.0020  | 0.5225 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0017  | 0.5842 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0057  | 0.5842 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/2   | 2        | 0.2364  | 0.7773 | 0.0      | 🔴      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation pme`.

| Source                               | Wall time | CPU                                      | Threads |
|:-------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260802_002454_pme_full.csv` | 31.3 min  | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
