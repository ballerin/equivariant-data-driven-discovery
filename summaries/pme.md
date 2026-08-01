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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best        | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0011     | multigraph/RFE     | 0.0019    | WSINDy/RFE      | 0/0 ✓     | 0.0010      | 0.0019     | 0.0022       |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0008     | multigraph/RFE     | 0.0051    | PolyJet/RFE     | 0/0 ✓     | 0.0008      | 0.0051     | 0.0052       |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0008     | multigraph/MSTLS   | 0.0006    | WSINDy/RFE      | 0/0 ✓     | 0.0007      | 0.0007     | 0.0010       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0007     | multigraph/RFE     | 0.0007    | WSINDy/RFE      | 0/0 ✓     | 0.0007      | 0.0007     | 0.0010       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0009     | multigraph/MSTLS   | 0.0010    | WSINDy/RFE      | 0/0 ✓     | 0.0008      | 0.0011     | 0.0014       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0009     | multigraph/STRidge | 0.0019    | PolyJet/STRidge | 0/0 ✓     | 0.0009      | 0.0019     | 0.0014       |
| 2%    | pointwise     | spectral   | —        | 21/0 ✗            | 0.9317     | —                  | 0.8527    | —               | 3/0 ✗     | 0.9254      | 0.8605     | 1.5402       |
| 2%    | pointwise     | polynomial | —        | 9/0 ✗             | 0.2217     | —                  | 0.8523    | —               | 3/0 ✗     | 0.1662      | 0.8545     | 0.8599       |
| 2%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0010     | multigraph/STRidge | 0.2406    | WSINDy/RFE      | 0/0 ✓     | 0.0016      | 0.2405     | 0.2406       |
| 2%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0009     | multigraph/STRidge | 0.2405    | WSINDy/RFE      | 0/0 ✓     | 0.0017      | 0.2405     | 0.2405       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0014     | multigraph/STRidge | 0.2795    | WSINDy/STRidge  | 0/0 ✓     | 0.0015      | 0.2795     | 0.2796       |
| 2%    | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0012     | multigraph/MSTLS   | 0.4409    | WSINDy/RFE      | 0/0 ✓     | 0.0010      | 0.4409     | 0.2794       |
| 5%    | pointwise     | spectral   | —        | 25/0 ✗            | 1.0014     | —                  | 0.9700    | —               | 36/0 ✗    | 1.0065      | 0.9701     | 1.7914       |
| 5%    | pointwise     | polynomial | —        | 4/0 ✗             | 0.5178     | —                  | 0.9716    | —               | 3/0 ✗     | 0.3398      | 0.9759     | 0.9825       |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0022     | multigraph/RFE     | 0.5229    | WSINDy/STRidge  | 0/0 ✓     | 0.0043      | 0.5225     | 0.5231       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0018     | multigraph/STRidge | 0.5225    | WSINDy/STRidge  | 0/0 ✓     | 0.0060      | 0.5225     | 0.5225       |
| 5%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0023     | multigraph/MSTLS   | 0.5848    | WSINDy/STRidge  | 0/0 ✓     | 0.0039      | 0.5842     | 0.5850       |
| 5%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0015     | multigraph/STRidge | 0.5842    | WSINDy/STRidge  | 0/0 ✓     | 0.0052      | 0.5842     | 0.5842       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best   |
|:----------------------|:--------------------------------|:------|--------:|---------:|-----------:|
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0007  | 0.1      | —          |
| WSINDy + STRidge      | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0007  | 0.1      | +2.8%      |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.0      | +3.2%      |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 14.4     | +3.2%      |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.2      | +3.2%      |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +3.2%      |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 0.1      | +6.3%      |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0007  | 2.8      | +6.3%      |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.8360  | 0.3      | +118725.8% |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 2/2 ✗ | 0.3963  | 0.2      | ✗          |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 2/2 ✗ | 1.1109  | 1.0      | ✗          |
| PDE-FIND + RFE        | spectral, pointwise             | 2/4 ✗ | 0.2585  | 0.4      | ✗          |

### 2% noise

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:----------------------|:--------------------------------|:------|--------:|---------:|----------:|
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0009  | 0.1      | —         |
| WSINDy + RFE          | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0010  | 0.0      | +0.4%     |
| Multigraphs + MSTLS   | polynomial, weak (IBP), m=8     | 0/0 ✓ | 0.0012  | 0.2      | +29.1%    |
| WSINDy + STRidge      | spectral, weak (IBP), m=12      | 0/0 ✓ | 0.0015  | 0.1      | +62.1%    |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0020  | 0.0      | +106.8%   |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8    | 0/0 ✓ | 0.0023  | 27.7     | +144.3%   |
| PolyJet + STRidge     | polynomial, weak (IBP), m=12    | 0/0 ✓ | 0.0084  | 2.7      | +786.7%   |
| WSINDy + MSTLS        | spectral, weak (no IBP), m=4    | 0/0 ✓ | 0.5144  | 0.7      | +54119.4% |
| PolyJet + RFE         | polynomial, weak (IBP), m=12    | 0/1 ✗ | 0.6061  | 0.1      | ✗         |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 1/2 ✗ | 0.3935  | 0.1      | ✗         |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=4    | 3/2 ✗ | 0.7039  | 1.3      | ✗         |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=4  | 2/4 ✗ | 0.3407  | 0.1      | ✗         |

### 5% noise

| Method                | Best conf.                     | fp/fn | E_coeff | Time (s) | vs. best |
|:----------------------|:-------------------------------|:------|--------:|---------:|---------:|
| Multigraphs + STRidge | polynomial, weak (IBP), m=12   | 0/0 ✓ | 0.0015  | 0.1      | —        |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12  | 0/0 ✓ | 0.0022  | 0.0      | +40.8%   |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12  | 0/0 ✓ | 0.0022  | 0.1      | +40.8%   |
| WSINDy + STRidge      | spectral, weak (IBP), m=12     | 0/0 ✓ | 0.0039  | 0.1      | +152.0%  |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=8   | 0/0 ✓ | 0.0045  | 25.9     | +193.8%  |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=4 | 0/0 ✓ | 0.0441  | 0.7      | +2767.4% |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=8 | 0/0 ✓ | 0.0651  | 3.0      | +4133.1% |
| WSINDy + RFE          | spectral, weak (IBP), m=8      | 0/2 ✗ | 0.4909  | 0.0      | ✗        |
| PDE-FIND + MSTLS      | spectral, weak (IBP), m=8      | 2/2 ✗ | 0.3327  | 1.0      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=8   | 2/2 ✗ | 0.5418  | 0.2      | ✗        |
| PolyJet + RFE         | spectral, weak (IBP), m=8      | 0/3 ✗ | 0.5828  | 0.4      | ✗        |
| PDE-FIND + RFE        | polynomial, pointwise          | 1/4 ✗ | 0.7899  | 0.3      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 7/2   | 14       | 0.7859  | 0.1747 | 0.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/4   | 2        | 0.2585  | 0.2565 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.001   | 19/2  | 26       | 1.2016  | 0.8578 | 0.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 2/4   | 2        | 0.9285  | 0.8606 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 23/2  | 37       | 1.1186  | 0.9691 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/4   | 1        | 0.9901  | 0.9854 | 0.2      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 0/0   | 2        | 0.0011  | 0.0019 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 2        | 0.0010  | 0.0019 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.01    | 3/0   | 5        | 0.9254  | 0.8605 | 0.2      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 2/4   | 2        | 0.9314  | 0.8647 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 36/0  | 22       | 1.0065  | 0.9701 | 0.2      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0011  | 0.0019 | 3.0      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 4        | 0.0011  | 0.0019 | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 5/0   | 9        | 0.9262  | 0.8556 | 5.0      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/4   | 2        | 0.9285  | 0.8606 | 0.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.001   | 75/0  | 79       | 1.0066  | 0.9666 | 8.5      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/4   | 1        | 0.9901  | 0.9854 | 0.4      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 0/0   | 2        | 0.0011  | 0.0019 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 2        | 0.0011  | 0.0019 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.01    | 21/0  | 11       | 0.9317  | 0.8527 | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 25/0  | 13       | 1.0014  | 0.9700 | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 2/2   | 6        | 2.2731  | 0.4066 | 0.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/4   | 2        | 0.2590  | 0.2565 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.01    | 2/2   | 6        | 2.3180  | 0.8796 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 2/4   | 2        | 0.3413  | 0.8606 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.001   | 32/2  | 45       | 1.9205  | 0.9700 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 1/4   | 1        | 0.7899  | 0.9854 | 0.3      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.001   | 0/0   | 2        | 0.0013  | 0.0057 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 2        | 0.0012  | 0.0057 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 3/0   | 5        | 0.1662  | 0.8545 | 0.2      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 0/2   | 1        | 0.5852  | 0.9065 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 3/0   | 4        | 0.3398  | 0.9759 | 0.2      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 1/4   | 1        | 0.7899  | 0.9854 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0009  | 0.0051 | 3.2      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 4        | 0.0008  | 0.0051 | 0.3      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.01    | 20/0  | 24       | 0.2833  | 0.8488 | 4.7      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 2/4   | 2        | 0.3413  | 0.8606 | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 4/0   | 8        | 0.5201  | 0.9716 | 3.7      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 1/4   | 1        | 0.7899  | 0.9854 | 0.3      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.01    | 0/0   | 2        | 0.0009  | 0.0051 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 2        | 0.0008  | 0.0051 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.01    | 9/0   | 7        | 0.2217  | 0.8523 | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 2/4   | 2        | 0.3413  | 0.8606 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.01    | 4/0   | 5        | 0.5178  | 0.9716 | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 1/4   | 1        | 0.7899  | 0.9854 | 0.0      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | no  | 2/2   | 4        | 0.3963  | 0.3115 | 0.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 2/2   | 8        | 1.1109  | 0.1749 | 1.0      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 2/4   | 2        | 0.2585  | 0.2518 | 0.1      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 12 | yes | 5/2   | 7        | 0.3008  | 0.1790 | 0.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 2/2   | 8        | 1.1137  | 0.1702 | 0.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/4   | 2        | 0.2586  | 0.2485 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 4  | no  | 3/2   | 7        | 0.7177  | 0.5868 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 4  | no  | 3/2   | 7        | 0.7039  | 0.5876 | 1.3      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 2/4   | 2        | 0.9286  | 0.5949 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 3/2   | 7        | 1.1563  | 0.4804 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 4/2   | 8        | 0.8932  | 0.4824 | 1.0      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.9287  | 0.4936 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 8  | no  | 2/2   | 6        | 0.5418  | 0.6662 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | no  | 2/2   | 6        | 0.6894  | 0.6660 | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/4   | 1        | 0.9902  | 0.9255 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 3/2   | 6        | 0.6217  | 0.6074 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 8  | yes | 2/2   | 6        | 0.3327  | 0.7773 | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.9805  | 0.7809 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.8369  | 0.8354 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0008  | 0.0011 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0017  | 0.2406 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.5144  | 0.7166 | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0016  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0015  | 0.2795 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.6468  | 0.7323 | 0.5      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0021  | 0.2795 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0043  | 0.5225 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 2        | 0.4008  | 0.7141 | 0.5      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 12 | no  | 0/2   | 1        | 0.4924  | 0.7166 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0039  | 0.5842 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.1729  | 0.7787 | 0.5      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 0/2   | 1        | 0.4909  | 0.8520 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 2.8      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 13.6     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0010  | 0.0010 | 2.8      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0009  | 0.0010 | 13.9     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0009  | 0.0010 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0310  | 0.3304 | 2.9      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0023  | 0.3303 | 27.7     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/1   | 3        | 1.2472  | 0.4159 | 0.5      | 🔴      |
| PolyJet     | 2%    | STRidge   | 0.1     | 8  | yes | 0/0   | 4        | 0.0559  | 0.4422 | 3.6      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0023  | 0.2795 | 19.6     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 8  | yes | 0/1   | 3        | 1.1067  | 0.5141 | 0.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.1019  | 0.5229 | 2.9      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0045  | 0.6475 | 25.9     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | no  | 2/3   | 3        | 25.3618 | 0.7893 | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 8  | yes | 1/0   | 5        | 0.3055  | 0.7683 | 3.6      | 🔴      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0086  | 0.5848 | 22.3     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/3   | 1        | 0.5828  | 0.8897 | 0.4      | 🔴      |
| Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0010  | 0.0010 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0009  | 0.0010 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0009  | 0.0010 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0010  | 0.2406 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0020  | 0.2406 | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0020  | 0.2406 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0014  | 0.2795 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0021  | 0.2795 | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 2        | 0.0021  | 0.2795 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 2        | 0.0032  | 0.6476 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0022  | 0.5229 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0022  | 0.5229 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0056  | 0.5848 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0023  | 0.5848 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 0/0   | 2        | 0.0250  | 0.7696 | 0.0      | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 1.0     | 8  | no  | 4/2   | 10       | 0.8816  | 0.1708 | 0.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 4  | no  | 2/2   | 8        | 1.1227  | 0.1769 | 1.1      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 2/4   | 2        | 0.2590  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 8  | yes | 3/2   | 9        | 1.7672  | 0.1667 | 0.2      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 8  | yes | 2/2   | 8        | 1.1260  | 0.1715 | 0.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/4   | 2        | 0.2591  | 0.2485 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 10.0    | 12 | no  | 1/2   | 3        | 0.3935  | 0.4377 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 4  | no  | 4/2   | 8        | 1.7944  | 0.5872 | 1.4      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 4  | no  | 2/4   | 2        | 0.3407  | 0.5949 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 8  | yes | 4/2   | 7        | 0.7023  | 0.4899 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 8  | yes | 4/2   | 8        | 1.8785  | 0.4824 | 1.1      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.3408  | 0.4936 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 4  | no  | 3/2   | 7        | 1.2938  | 0.8568 | 0.3      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 4  | no  | 2/2   | 6        | 1.0094  | 0.8569 | 1.4      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 1/4   | 1        | 0.7902  | 0.9255 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 8  | yes | 3/2   | 8        | 1.3893  | 0.7772 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | MSTLS     | 0.75    | 8  | yes | 2/2   | 6        | 1.1156  | 0.7773 | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 2/4   | 2        | 0.5512  | 0.7809 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.01    | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 2        | 0.0009  | 0.0012 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.8369  | 0.8353 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0022  | 0.2406 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.5209  | 0.7119 | 0.7      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0017  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 2        | 0.0020  | 0.2795 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.6483  | 0.7311 | 0.5      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 8  | yes | 0/0   | 2        | 0.0010  | 0.4409 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0060  | 0.5225 | 0.2      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 4  | no  | 0/0   | 2        | 0.0441  | 0.8533 | 0.7      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 0/2   | 1        | 0.6069  | 0.9042 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0052  | 0.5842 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.1780  | 0.7766 | 0.6      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 0/2   | 1        | 0.6078  | 0.8513 | 0.1      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 2.8      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 14.4     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 8  | yes | 0/0   | 4        | 0.0009  | 0.0019 | 3.0      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0009  | 0.0011 | 14.5     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0009  | 0.0011 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0151  | 0.3303 | 3.0      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0033  | 0.3302 | 27.1     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 8  | no  | 0/2   | 2        | 0.3903  | 0.4647 | 0.3      | 🔴      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0084  | 0.2795 | 2.7      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0033  | 0.2794 | 21.8     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/1   | 3        | 0.6061  | 0.3767 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0651  | 0.6471 | 3.0      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | no  | 0/0   | 4        | 0.0078  | 0.6469 | 25.5     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | no  | 0/3   | 1        | 0.6747  | 0.8245 | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.1214  | 0.7671 | 3.5      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 4        | 0.0272  | 0.7668 | 22.9     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/3   | 1        | 0.6761  | 0.8780 | 0.3      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 0.01    | 8  | yes | 0/0   | 2        | 0.0009  | 0.0019 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 2        | 0.0009  | 0.0011 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 2        | 0.0009  | 0.2405 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0028  | 0.2405 | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0028  | 0.2405 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 8  | yes | 0/0   | 2        | 0.0017  | 0.4409 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 8  | yes | 0/0   | 2        | 0.0012  | 0.4409 | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 8  | yes | 3/4   | 1        | 0.6784  | 0.7657 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 2        | 0.0018  | 0.5225 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 12 | no  | 0/0   | 2        | 0.0056  | 0.5225 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 2        | 0.0015  | 0.5842 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 2        | 0.0056  | 0.5842 | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 8  | yes | 1/2   | 2        | 0.5390  | 0.7773 | 0.0      | 🔴      |


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
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0019       | 1215766.5 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0014       | 1837081.4 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0021       | 610506.4  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0014       | 1106284.3 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.5638       | 5461.2    |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.3303       | 14793.5   |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.2406       | 25197.9   |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.5597       | 26669.9   |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3302       | 75261.8   |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.2405       | 169466.9  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.4421       | 7367.5    |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.2796       | 14498.1   |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.4409       | 31979.0   |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.2794       | 78666.4   |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.8612       | 2218.6    |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.6476       | 5152.9    |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.5231       | 11056.3   |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.8536       | 9035.7    |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.6469       | 30598.9   |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.5225       | 63503.4   |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.7698       | 2723.7    |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.5850       | 5539.1    |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.7669       | 12000.1   |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.5842       | 31477.2   |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation pme`.

| Source                               | Wall time | CPU                                      | Threads |
|:-------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260730_030731_pme_full.csv` | 104.3 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | not set |
