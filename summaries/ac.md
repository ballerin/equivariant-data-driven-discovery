# Allen-Cahn equation

The equation under consideration is

$$f_t = f - f^3 + \varepsilon^2\Delta f, \qquad \varepsilon = 0.15$$

which can be written in terms of $O(2)$-equivariant terms

| Operator   | True coefficient           |
|:-----------|:---------------------------|
| $f$        | $+1.0$                     |
| $f^3$      | $-1.0$                     |
| $\Delta f$ | $+\varepsilon^2 = +0.0225$ |

or standard terms (PolyJet basis)

| Operator          | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f$               |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^3$             |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |

Allen-Cahn contains only pure monomials ($f$, $f^3$) and pure second derivatives ($f_{xx}$, $f_{yy}$).  All four libraries can represent every term, and this is considered the simplest equation in the benchmark set.

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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best         | Mg E_data | Std best         | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-------------------|:----------|:-----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0080     | multigraph/RFE     | 0.0664    | PDE-FIND/RFE     | 0/0 ✓     | 0.0080      | 0.0664     | 0.0669       |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0113     | multigraph/RFE     | 0.1133    | PDE-FIND/RFE     | 0/0 ✓     | 0.0114      | 0.1133     | 0.1139       |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge | 0.0032    | PolyJet/STRidge  | 0/0 ✓     | 0.0048      | 0.0032     | 0.0057       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge | 0.0032    | PDE-FIND/STRidge | 0/0 ✓     | 0.0048      | 0.0032     | 0.0058       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0050     | multigraph/STRidge | 0.0034    | PDE-FIND/STRidge | 0/0 ✓     | 0.0050      | 0.0034     | 0.0061       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0050     | multigraph/STRidge | 0.0034    | PDE-FIND/STRidge | 0/0 ✓     | 0.0050      | 0.0034     | 0.0061       |
| 2%    | pointwise     | spectral   | —        | 3/0 ✗             | 0.6957     | —                  | 0.3221    | —                | 2/0 ✗     | 0.6989      | 0.3047     | 1.0395       |
| 2%    | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0272     | multigraph/STRidge | 0.2545    | PolyJet/STRidge  | 0/0 ✓     | 0.0280      | 0.2545     | 0.2560       |
| 2%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0050     | multigraph/STRidge | 0.0217    | PolyJet/STRidge  | 0/0 ✓     | 0.0050      | 0.0217     | 0.0222       |
| 2%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/RFE     | 0.0217    | PDE-FIND/RFE     | 0/0 ✓     | 0.0048      | 0.0217     | 0.0222       |
| 2%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0054     | multigraph/STRidge | 0.0287    | PDE-FIND/STRidge | 0/0 ✓     | 0.0054      | 0.0287     | 0.0292       |
| 2%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0051     | multigraph/RFE     | 0.0287    | WSINDy/RFE       | 0/0 ✓     | 0.0051      | 0.0287     | 0.0292       |
| 5%    | pointwise     | spectral   | —        | 3/0 ✗             | 0.9248     | —                  | 0.5486    | —                | 2/0 ✗     | 0.9250      | 0.5437     | 2.3614       |
| 5%    | pointwise     | polynomial | 🟢        | 0/0 ✓             | 0.1447     | multigraph/STRidge | 0.5150    | —                | 4/0 ✗     | 0.1488      | 0.5140     | 0.5325       |
| 5%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0070     | multigraph/STRidge | 0.0536    | PolyJet/STRidge  | 0/0 ✓     | 0.0069      | 0.0536     | 0.0540       |
| 5%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0065     | multigraph/MSTLS   | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0064      | 0.0536     | 0.0540       |
| 5%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0088     | multigraph/STRidge | 0.0713    | PolyJet/STRidge  | 0/0 ✓     | 0.0088      | 0.0713     | 0.0717       |
| 5%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0071     | multigraph/MSTLS   | 0.0713    | PDE-FIND/MSTLS   | 0/0 ✓     | 0.0071      | 0.0713     | 0.0717       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:----------------------|:--------------------------------|:------|--------:|---------:|----------:|
| Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | —         |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 2.8      | +0.4%     |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.3      | +0.4%     |
| WSINDy + STRidge      | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.4%     |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.2      | +0.8%     |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.8%     |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 15.1     | +1.1%     |
| PDE-FIND + MSTLS      | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 2.3      | +1.1%     |
| PolyJet + RFE         | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +1.1%     |
| PDE-FIND + RFE        | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +1.1%     |
| WSINDy + RFE          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +1.1%     |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 2.8299  | 0.3      | +59028.8% |

### 2% noise

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:----------------------|:--------------------------------|:------|--------:|---------:|----------:|
| Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | —         |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PolyJet + RFE         | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.1      | +0.2%     |
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PDE-FIND + MSTLS      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 0.7      | +0.2%     |
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0048  | 15.6     | +0.2%     |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0049  | 0.1      | +2.0%     |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0049  | 2.8      | +2.4%     |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0049  | 0.2      | +2.4%     |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0049  | 0.1      | +2.4%     |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 2.8777  | 0.3      | +59464.8% |

### 5% noise

| Method                | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:----------------------|:--------------------------------|:------|--------:|---------:|----------:|
| PolyJet + MSTLS       | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0064  | 19.0     | —         |
| PDE-FIND + MSTLS      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0064  | 0.6      | +0.0%     |
| Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.1      | +0.3%     |
| PolyJet + STRidge     | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 2.9      | +2.9%     |
| PDE-FIND + STRidge    | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.2      | +2.9%     |
| WSINDy + STRidge      | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.1      | +2.9%     |
| Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.1      | +3.1%     |
| Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0068  | 0.0      | +5.4%     |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 3.1029  | 0.3      | +48030.2% |
| PDE-FIND + RFE        | polynomial, weak (no IBP), m=4  | 0/2 ✗ | 0.3967  | 0.1      | ✗         |
| WSINDy + RFE          | polynomial, weak (no IBP), m=4  | 0/2 ✗ | 0.3967  | 0.0      | ✗         |
| PolyJet + RFE         | polynomial, weak (no IBP), m=4  | 0/2 ✗ | 0.3967  | 0.3      | ✗         |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0084  | 0.0664 | 0.2      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.3      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 2/0   | 6        | 0.6989  | 0.3047 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | RFE       |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.1     | 2/0   | 6        | 0.9250  | 0.5437 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.2      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0084  | 0.0664 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 2/0   | 6        | 0.6989  | 0.3047 | 0.1      | 🔴      |
| WSINDy      | 2%    | RFE       |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 2/0   | 6        | 0.9250  | 0.5437 | 0.2      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0084  | 0.0664 | 3.9      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 2/0   | 6        | 0.6989  | 0.3047 | 3.5      | 🔴      |
| PolyJet     | 2%    | RFE       |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 17/0  | 21       | 0.9271  | 0.5112 | 4.8      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.4      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0084  | 0.0664 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 3/0   | 4        | 0.6957  | 0.3221 | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE       |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 3/0   | 4        | 0.9248  | 0.5486 | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE       |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0116  | 0.1133 | 0.2      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 0/0   | 4        | 0.0114  | 0.1133 | 0.2      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 0/0   | 4        | 0.0280  | 0.2545 | 0.2      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 0/2   | 2        | 0.3414  | 0.4072 | 0.4      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 0.01    | 4/0   | 8        | 0.1488  | 0.5140 | 0.3      | 🔴      |
| PDE-FIND    | 5%    | RFE       |         | 0/2   | 2        | 0.3973  | 0.5669 | 0.4      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 0/0   | 4        | 0.0116  | 0.1133 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 4        | 0.0114  | 0.1133 | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 0/0   | 4        | 0.0280  | 0.2545 | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 0/2   | 2        | 0.3414  | 0.4072 | 0.1      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.01    | 8/0   | 9        | 0.1780  | 0.5085 | 0.2      | 🔴      |
| WSINDy      | 5%    | RFE       |         | 0/2   | 2        | 0.3973  | 0.5669 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0116  | 0.1133 | 3.2      | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 4        | 0.0114  | 0.1133 | 0.4      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 0/0   | 4        | 0.0280  | 0.2545 | 3.1      | 🟢      |
| PolyJet     | 2%    | RFE       |         | 0/2   | 2        | 0.3414  | 0.4072 | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.01    | 12/0  | 16       | 0.2196  | 0.4925 | 4.9      | 🔴      |
| PolyJet     | 5%    | RFE       |         | 0/2   | 2        | 0.3973  | 0.5669 | 0.3      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0115  | 0.1133 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0113  | 0.1133 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 0/0   | 3        | 0.0272  | 0.2545 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 0/2   | 2        | 0.3414  | 0.4072 | 0.0      | 🔴      |
| Multigraphs | 5%    | STRidge   | 0.1     | 0/0   | 3        | 0.1447  | 0.5150 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 0/2   | 2        | 0.3973  | 0.5669 | 0.0      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 2.3      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.6      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0050  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.6      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0054  | 0.0287 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0055  | 0.0287 | 0.5      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0055  | 0.0287 | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0078  | 0.0536 | 0.6      | 🟢      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0088  | 0.0713 | 0.2      | 🟢      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0091  | 0.0713 | 0.6      | 🟢      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.9551  | 0.3001 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8796  | 2.8813 | 0.4      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0766  | 2.7381 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 4.9910  | 2.8793 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0054  | 0.0287 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 6.0202  | 2.7348 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0055  | 0.0287 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 6.0213  | 2.8685 | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0088  | 0.0713 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 7.3772  | 2.7195 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.9551  | 0.3001 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 2.8      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 15.1     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 3.2      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 13.6     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0050  | 0.0217 | 3.5      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 15.8     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0054  | 0.0287 | 3.1      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0055  | 0.0287 | 14.7     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0055  | 0.0287 | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0069  | 0.0536 | 2.9      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0078  | 0.0536 | 16.9     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.4      | 🔴      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0088  | 0.0713 | 2.6      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0091  | 0.0713 | 17.1     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.9551  | 0.3001 | 0.3      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0054  | 0.0287 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0055  | 0.0287 | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0055  | 0.0287 | 0.1      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0070  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0078  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0103  | 0.0770 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0088  | 0.0713 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0091  | 0.0713 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0091  | 0.0713 | 0.0      | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 1.5      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND    | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.2      | 🟢      |
| PDE-FIND    | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.5      | 🟢      |
| PDE-FIND    | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.7      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0053  | 0.0287 | 0.2      | 🟢      |
| PDE-FIND    | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.5      | 🟢      |
| PDE-FIND    | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0064  | 0.0536 | 0.6      | 🟢      |
| PDE-FIND    | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3967  | 0.3345 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0076  | 0.0713 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0071  | 0.0713 | 0.6      | 🟢      |
| PDE-FIND    | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3997  | 0.3001 | 0.1      | 🔴      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8299  | 2.8812 | 0.3      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0014  | 2.7381 | 0.4      | 🟢      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8777  | 2.8792 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0053  | 0.0287 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0742  | 2.7348 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 3.1029  | 2.8685 | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3967  | 0.3345 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0076  | 0.0713 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.4085  | 2.7195 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3997  | 0.3001 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 3.3      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 14.4     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 2.9      | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 13.5     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 2.8      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 15.6     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0053  | 0.0287 | 3.0      | 🟢      |
| PolyJet     | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 14.4     | 🟢      |
| PolyJet     | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 2.9      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0064  | 0.0536 | 19.0     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3967  | 0.3345 | 0.3      | 🔴      |
| PolyJet     | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0076  | 0.0713 | 4.4      | 🟢      |
| PolyJet     | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0071  | 0.0713 | 15.7     | 🟢      |
| PolyJet     | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3997  | 0.3001 | 0.3      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.2      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0053  | 0.0287 | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.2      | 🟢      |
| Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0065  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0068  | 0.0770 | 0.0      | 🟢      |
| Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0076  | 0.0713 | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0071  | 0.0713 | 0.2      | 🟢      |
| Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0071  | 0.0713 | 0.0      | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0070       | 1742.0  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0061       | 1475.2  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0057       | 2683.0  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0168       | 1704.3  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0063       | 1431.7  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0058       | 2507.3  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0066       | 1378.1  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0061       | 1745.6  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0066       | 1396.1  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0061       | 1678.4  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.0713       | 2330.4  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.0315       | 2821.4  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.0222       | 4037.5  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.0655       | 1731.0  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.0315       | 1425.0  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.0222       | 2331.3  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.0495       | 2402.4  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.0292       | 3028.3  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.0495       | 1423.1  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.0292       | 1665.4  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.1759       | 10382.8 |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.0773       | 13252.1 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.0540       | 19563.0 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.1578       | 2043.8  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.0773       | 1387.5  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.0540       | 1856.6  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.1224       | 10811.9 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.0717       | 14181.2 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.1224       | 1694.1  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.0717       | 1586.1  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ac`.

| Source                              | Wall time | CPU                                      | Threads |
|:------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260729_180931_ac_full.csv` | 114.3 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | not set |
