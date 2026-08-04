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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best       | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:---------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0010     | mod_multigraph/RFE     | 0.0019    | WSINDy/RFE     | 0/0 ✓     | 0.0010      | 0.0019     | 0.0022       |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0008     | mod_multigraph/RFE     | 0.0051    | PolyJet/RFE    | 0/0 ✓     | 0.0008      | 0.0051     | 0.0052       |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0007     | mod_multigraph/MSTLS   | 0.0007    | WSINDy/RFE     | 0/0 ✓     | 0.0007      | 0.0007     | 0.0010       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0007     | mod_multigraph/RFE     | 0.0007    | WSINDy/RFE     | 0/0 ✓     | 0.0007      | 0.0007     | 0.0010       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0007     | mod_multigraph/MSTLS   | 0.0007    | WSINDy/RFE     | 0/0 ✓     | 0.0007      | 0.0007     | 0.0010       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0007     | mod_multigraph/MSTLS   | 0.0007    | WSINDy/RFE     | 0/0 ✓     | 0.0007      | 0.0007     | 0.0010       |
| 2%    | pointwise     | spectral   | —        | 10/0 ✗            | 0.2595     | —                      | 0.8545    | —              | 12/0 ✗    | 0.2552      | 0.8601     | 1.5402       |
| 2%    | pointwise     | polynomial | —        | 5/0 ✗             | 0.0778     | —                      | 0.8587    | —              | 3/0 ✗     | 0.0388      | 0.8533     | 0.8599       |
| 2%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0006     | mod_multigraph/RFE     | 0.2406    | WSINDy/RFE     | 0/0 ✓     | 0.0017      | 0.2405     | 0.2406       |
| 2%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0005     | mod_multigraph/STRidge | 0.2405    | WSINDy/RFE     | 0/0 ✓     | 0.0017      | 0.2405     | 0.2405       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0000     | mod_multigraph/STRidge | 0.2405    | WSINDy/RFE     | 0/0 ✓     | 0.0017      | 0.2405     | 0.2406       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0000     | mod_multigraph/STRidge | 0.2405    | WSINDy/RFE     | 0/0 ✓     | 0.0017      | 0.2405     | 0.2405       |
| 5%    | pointwise     | spectral   | —        | 32/0 ✗            | 0.7804     | —                      | 0.9699    | —              | 22/0 ✗    | 0.8836      | 0.9702     | 1.7914       |
| 5%    | pointwise     | polynomial | —        | 4/0 ✗             | 0.1883     | —                      | 0.9716    | —              | 7/0 ✗     | 0.2165      | 0.9722     | 0.9825       |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0026     | mod_multigraph/MSTLS   | 0.6476    | WSINDy/STRidge | 0/0 ✓     | 0.0056      | 0.5225     | 0.5231       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0017     | mod_multigraph/STRidge | 0.5225    | WSINDy/STRidge | 0/0 ✓     | 0.0056      | 0.5225     | 0.5225       |
| 5%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0013     | mod_multigraph/STRidge | 0.5225    | WSINDy/STRidge | 0/0 ✓     | 0.0056      | 0.5225     | 0.5231       |
| 5%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0013     | mod_multigraph/STRidge | 0.5225    | WSINDy/STRidge | 0/0 ✓     | 0.0056      | 0.5225     | 0.5225       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best   |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|-----------:|
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.0      | —          |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +0.0%      |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0007  | 0.0      | +6.2%      |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0007  | 0.1      | +8.8%      |
| WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0007  | 0.2      | +9.1%      |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.0      | +9.6%      |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.8      | +9.6%      |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +9.6%      |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +9.6%      |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +12.8%     |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.5      | +12.8%     |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.8360  | 0.2      | +126088.9% |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=12   | 2/2 ✗ | 0.3963  | 0.3      | ✗          |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4    | 6/2 ✗ | 0.1893  | 0.8      | ✗          |
| PDE-FIND + RFE                 | spectral, pointwise             | 2/4 ✗ | 0.2585  | 0.1      | ✗          |

### 2% noise

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best    |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|------------:|
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0000  | 0.1      | —           |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0006  | 0.0      | +2549.3%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0006  | 0.1      | +2549.3%    |
| Multigraphs + STRidge          | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0010  | 0.1      | +4275.6%    |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0017  | 0.0      | +7464.8%    |
| WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0022  | 0.2      | +9901.4%    |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0028  | 0.0      | +12850.9%   |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0028  | 0.1      | +12850.9%   |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0034  | 1.6      | +15157.7%   |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0068  | 0.5      | +30992.2%   |
| WSINDy + MSTLS                 | spectral, weak (no IBP), m=4    | 0/0 ✓ | 0.5144  | 0.5      | +2341710.0% |
| PolyJet + RFE                  | polynomial, weak (IBP), m=8     | 0/1 ✗ | 0.2975  | 0.1      | ✗           |
| PDE-FIND + STRidge             | polynomial, weak (no IBP), m=12 | 1/2 ✗ | 0.4002  | 0.2      | ✗           |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=4    | 5/2 ✗ | 0.2326  | 0.9      | ✗           |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=4    | 2/4 ✗ | 0.2586  | 0.0      | ✗           |

### 5% noise

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|---------:|
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0013  | 0.1      | —        |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0020  | 0.1      | +55.2%   |
| Modified Multigraphs + MSTLS   | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0025  | 0.1      | +93.0%   |
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0025  | 0.0      | +93.0%   |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0056  | 0.1      | +334.6%  |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0056  | 0.0      | +334.6%  |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0056  | 0.2      | +336.0%  |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0059  | 1.6      | +364.3%  |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0193  | 0.6      | +1406.5% |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=4  | 0/0 ✓ | 0.0441  | 0.5      | +3346.2% |
| WSINDy + RFE                   | polynomial, weak (no IBP), m=8  | 0/2 ✗ | 0.5773  | 0.0      | ✗        |
| PDE-FIND + STRidge             | polynomial, weak (no IBP), m=4  | 2/2 ✗ | 0.5205  | 0.4      | ✗        |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=4  | 4/2 ✗ | 0.2377  | 0.9      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/3 ✗ | 0.6700  | 0.1      | ✗        |
| PDE-FIND + RFE                 | spectral, pointwise             | 1/4 ✗ | 0.7000  | 0.0      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 10/2  | 12       | 0.6329  | 0.2109 | 0.5      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 2/4   | 2        | 0.2585  | 0.2565 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.001   | 24/2  | 26       | 0.4869  | 0.8578 | 0.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 24/2  | 26       | 1.0370  | 0.9698 | 0.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 0/0   | 2        | 0.0011  | 0.0019 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 2        | 0.0010  | 0.0019 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.001   | 12/0  | 8        | 0.2552  | 0.8601 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 2/4   | 2        | 0.3045  | 0.8647 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 22/0  | 15       | 0.8836  | 0.9702 | 0.5      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 6/0   | 10       | 0.0016  | 0.0016 | 0.7      | 🔴      |
| PolyJet          | Clean | RFE       |         | 0/0   | 4        | 0.0011  | 0.0019 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 29/0  | 33       | 0.2409  | 0.8508 | 2.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 73/0  | 77       | 0.9578  | 0.9662 | 3.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 0/0   | 2        | 0.0011  | 0.0019 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 2        | 0.0011  | 0.0019 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 10/0  | 8        | 0.2595  | 0.8545 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.1      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 32/0  | 17       | 0.7804  | 0.9699 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 0/0   | 1        | 0.0011  | 0.0019 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 1        | 0.0010  | 0.0019 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 15/0  | 6        | 0.1855  | 0.8555 | 0.1      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 2/4   | 2        | 0.2595  | 0.8606 | 0.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 35/0  | 18       | 0.8047  | 0.9699 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/4   | 1        | 0.7000  | 0.9854 | 0.0      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.001   | 16/2  | 18       | 0.6129  | 0.1748 | 0.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 2/4   | 2        | 0.2590  | 0.2565 | 0.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 5/2   | 7        | 0.2541  | 0.8598 | 0.4      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.001   | 27/2  | 29       | 0.7600  | 0.9706 | 0.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 0/0   | 2        | 0.0013  | 0.0057 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 2        | 0.0012  | 0.0057 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 3/0   | 3        | 0.0388  | 0.8533 | 0.3      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 0/2   | 1        | 0.5784  | 0.9065 | 0.1      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.001   | 7/0   | 7        | 0.2165  | 0.9722 | 0.5      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 13/0  | 17       | 0.0029  | 0.0043 | 0.9      | 🔴      |
| PolyJet          | Clean | RFE       |         | 0/0   | 4        | 0.0008  | 0.0051 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 20/0  | 24       | 0.1182  | 0.8488 | 1.0      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.001   | 37/0  | 41       | 0.2905  | 0.9686 | 1.7      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 0/0   | 2        | 0.0009  | 0.0051 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 2        | 0.0008  | 0.0051 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 13/0  | 8        | 0.1117  | 0.8506 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.0      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 4/0   | 5        | 0.1883  | 0.9716 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.001   | 0/0   | 1        | 0.0009  | 0.0051 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 1        | 0.0008  | 0.0051 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 5/0   | 3        | 0.0778  | 0.8587 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 2/4   | 2        | 0.2600  | 0.8606 | 0.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 5/0   | 3        | 0.3166  | 0.9766 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/4   | 1        | 0.7001  | 0.9854 | 0.0      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 1.0     | 12 | no  | 2/2   | 4        | 0.3963  | 0.3115 | 0.3      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 6/2   | 8        | 0.1893  | 0.1749 | 0.8      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 2/4   | 2        | 0.2585  | 0.2518 | 0.1      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 1.0     | 12 | yes | 2/2   | 4        | 0.3963  | 0.3115 | 0.3      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 8  | yes | 6/2   | 8        | 0.2051  | 0.1604 | 0.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 2/4   | 2        | 0.2593  | 0.2445 | 0.1      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | no  | 3/2   | 5        | 0.5812  | 0.5944 | 0.3      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 4  | no  | 5/2   | 7        | 0.2326  | 0.5876 | 0.9      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 2/4   | 2        | 0.2586  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | yes | 3/2   | 5        | 0.5812  | 0.5944 | 0.3      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.75    | 8  | yes | 26/2  | 28       | 0.3100  | 0.3668 | 0.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2593  | 0.4044 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 8  | no  | 4/2   | 6        | 0.2576  | 0.6662 | 0.3      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | no  | 4/2   | 6        | 0.2617  | 0.6660 | 0.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/4   | 1        | 0.7004  | 0.9255 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 8  | yes | 4/2   | 6        | 0.2576  | 0.6662 | 0.3      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 8  | yes | 4/2   | 6        | 0.2617  | 0.6660 | 0.7      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2593  | 0.6743 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0022  | 0.2406 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.5144  | 0.7166 | 0.5      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0017  | 0.2405 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0022  | 0.2406 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.7339  | 0.7642 | 0.4      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0017  | 0.2405 | 0.0      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 2        | 0.4008  | 0.7141 | 0.4      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | no  | 0/2   | 1        | 0.5774  | 0.7821 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0056  | 0.5225 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.4008  | 0.7141 | 0.4      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 0/2   | 1        | 0.5774  | 0.7821 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0008  | 0.0006 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0008  | 0.0006 | 0.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0081  | 0.3304 | 0.7      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0042  | 0.3303 | 1.9      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/1   | 3        | 0.2979  | 0.4159 | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.0081  | 0.3304 | 0.7      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0042  | 0.3303 | 1.9      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/1   | 3        | 0.2979  | 0.4159 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | no  | 6/0   | 10       | 0.0564  | 0.5228 | 0.5      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0109  | 0.6475 | 1.4      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | no  | 2/3   | 3        | 0.6558  | 0.7893 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | yes | 6/0   | 10       | 0.0564  | 0.5228 | 0.5      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0109  | 0.6475 | 1.4      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 2/3   | 3        | 0.6558  | 0.7893 | 0.0      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 8  | no  | 0/0   | 2        | 0.0009  | 0.0011 | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs      | Clean | STRidge   | 0.01    | 8  | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.3      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0008  | 0.0006 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0026  | 0.2406 | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0038  | 0.2406 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0038  | 0.2406 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0026  | 0.2406 | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0038  | 0.2406 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0038  | 0.2406 | 0.1      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 2        | 0.0052  | 0.6476 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0107  | 0.5229 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0107  | 0.5229 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 8  | yes | 0/0   | 2        | 0.0052  | 0.6476 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0107  | 0.5229 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0107  | 0.5229 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 12 | no  | 0/0   | 1        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 1        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 1        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 1        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 1        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 1        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 1        | 0.0009  | 0.2406 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 1        | 0.0006  | 0.2406 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 1        | 0.0006  | 0.2406 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 1        | 0.0000  | 0.2405 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 1        | 0.0015  | 0.2405 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 1        | 0.0015  | 0.2405 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 1        | 0.0038  | 0.6476 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 1        | 0.0026  | 0.6476 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 1        | 0.0026  | 0.6476 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 1        | 0.0013  | 0.5225 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 1        | 0.0025  | 0.5225 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 1        | 0.0025  | 0.5225 | 0.0      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 4  | no  | 3/2   | 5        | 0.3673  | 0.3697 | 0.4      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 4  | no  | 6/2   | 8        | 0.1951  | 0.1769 | 0.6      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 4  | no  | 2/4   | 2        | 0.2590  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND         | Clean | STRidge   | 0.1     | 8  | yes | 2/2   | 4        | 0.5197  | 0.5382 | 0.3      | 🔴      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 8  | yes | 6/2   | 8        | 0.2147  | 0.1610 | 0.5      | 🔴      |
| PDE-FIND         | Clean | RFE       |         | 8  | yes | 2/4   | 2        | 0.2598  | 0.2445 | 0.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 10.0    | 12 | no  | 1/2   | 3        | 0.4002  | 0.4377 | 0.2      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 4  | no  | 6/2   | 8        | 0.2308  | 0.5872 | 0.8      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 4  | no  | 2/4   | 2        | 0.2590  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND         | 2%    | STRidge   | 10.0    | 12 | yes | 1/2   | 3        | 0.4002  | 0.4377 | 0.2      | 🔴      |
| PDE-FIND         | 2%    | MSTLS     | 0.5     | 8  | yes | 14/2  | 16       | 0.2367  | 0.3605 | 0.6      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2598  | 0.4044 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 4  | no  | 2/2   | 4        | 0.5205  | 0.8937 | 0.4      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 4  | no  | 4/2   | 6        | 0.2377  | 0.8569 | 0.9      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 1/4   | 1        | 0.7005  | 0.9255 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 10.0    | 12 | yes | 5/2   | 7        | 0.3393  | 0.5365 | 0.2      | 🔴      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 2/3   | 3        | 0.6679  | 0.6790 | 0.5      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.2598  | 0.6743 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.001   | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.001   | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0022  | 0.2406 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.5209  | 0.7119 | 0.5      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0017  | 0.2405 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0022  | 0.2406 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.7339  | 0.7641 | 0.4      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0017  | 0.2405 | 0.0      | 🟢      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.0441  | 0.8533 | 0.5      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | no  | 0/2   | 1        | 0.5773  | 0.7821 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0056  | 0.5225 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.4008  | 0.7140 | 0.4      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 0/2   | 1        | 0.5773  | 0.7821 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.8      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0007  | 0.0007 | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0007  | 0.0007 | 0.8      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0068  | 0.2406 | 0.5      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0034  | 0.3302 | 1.6      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/2   | 2        | 0.3503  | 0.4647 | 0.1      | 🔴      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0068  | 0.2406 | 0.5      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0034  | 0.3302 | 1.5      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/1   | 3        | 0.2975  | 0.4153 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0193  | 0.6471 | 0.6      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0059  | 0.6469 | 1.6      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | no  | 0/3   | 1        | 0.6700  | 0.8245 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.0193  | 0.6471 | 0.6      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0059  | 0.6469 | 1.7      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/3   | 1        | 0.6700  | 0.8245 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 2        | 0.0021  | 0.3302 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0028  | 0.2405 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0028  | 0.2405 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | yes | 0/0   | 2        | 0.0010  | 0.2405 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0028  | 0.2405 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0028  | 0.2405 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 2        | 0.0020  | 0.5225 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 2        | 0.0020  | 0.5225 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0056  | 0.5225 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0056  | 0.5225 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 8  | no  | 0/0   | 1        | 0.0008  | 0.0012 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 1        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 1        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 1        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 1        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 1        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 1        | 0.0005  | 0.2405 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 1        | 0.0020  | 0.2405 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 1        | 0.0020  | 0.2405 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 1        | 0.0000  | 0.2405 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 1        | 0.0015  | 0.2405 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 1        | 0.0015  | 0.2405 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 12 | no  | 0/0   | 1        | 0.0017  | 0.5225 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 1        | 0.0054  | 0.6469 | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 1        | 0.0054  | 0.6469 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 1        | 0.0013  | 0.5225 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 1        | 0.0025  | 0.5225 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 0/0   | 1        | 0.0031  | 0.6470 | 0.0      | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)   |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:----------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0020       | 1144526.1 |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0014       | 1933914.7 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0010       | 3063426.0 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0026       | 504516.0  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0015       | 1205386.4 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0010       | 1995810.1 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0014       | 1940044.0 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0010       | 3065203.8 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0015       | 1156385.6 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0010       | 1983785.4 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.5638       | 5461.2    |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.3303       | 14793.5   |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.2406       | 25197.9   |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.5597       | 26669.9   |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3302       | 75261.8   |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.2405       | 169466.9  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.3303       | 14794.1   |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.2406       | 25197.6   |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.3302       | 73762.5   |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.2405       | 169303.6  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.8612       | 2218.6    |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.6476       | 5152.9    |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.5231       | 11056.3   |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.8536       | 9035.7    |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.6469       | 30598.9   |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.5225       | 63503.4   |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.6476       | 5133.2    |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.5231       | 11050.3   |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.6469       | 30382.2   |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.5225       | 63501.7   |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation pme`.

| Source                                   | Wall time | CPU | Threads |
|:-----------------------------------------|----------:|:----|:-------:|
| `sweep_no_ac_ch_pme_20260803_170832.csv` | ? min     | ?   | ?       |
