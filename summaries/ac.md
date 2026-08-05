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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best         | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:-----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0080     | multigraph/RFE         | 0.0664    | WSINDy/RFE       | 0/0 ✓     | 0.0080      | 0.0664     | 0.0669       |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0114     | mod_multigraph/RFE     | 0.1133    | WSINDy/RFE       | 0/0 ✓     | 0.0115      | 0.1133     | 0.1139       |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge     | 0.0032    | PolyJet/STRidge  | 0/0 ✓     | 0.0048      | 0.0032     | 0.0057       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048     | mod_multigraph/STRidge | 0.0032    | PolyJet/STRidge  | 0/0 ✓     | 0.0048      | 0.0032     | 0.0058       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge     | 0.0032    | WSINDy/STRidge   | 0/0 ✓     | 0.0048      | 0.0032     | 0.0057       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge     | 0.0032    | PDE-FIND/STRidge | 0/0 ✓     | 0.0048      | 0.0032     | 0.0057       |
| 2%    | pointwise     | spectral   | —        | 40/0 ✗            | 0.1497     | —                      | 0.2885    | —                | 2/0 ✗     | 0.1340      | 0.3047     | 1.0395       |
| 2%    | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0269     | multigraph/STRidge     | 0.2545    | PDE-FIND/STRidge | 0/0 ✓     | 0.0278      | 0.2545     | 0.2560       |
| 2%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0048     | mod_multigraph/RFE     | 0.0217    | WSINDy/RFE       | 0/0 ✓     | 0.0048      | 0.0217     | 0.0222       |
| 2%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/RFE         | 0.0217    | WSINDy/RFE       | 0/0 ✓     | 0.0048      | 0.0217     | 0.0222       |
| 2%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0048     | mod_multigraph/RFE     | 0.0217    | WSINDy/RFE       | 0/0 ✓     | 0.0048      | 0.0217     | 0.0222       |
| 2%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/RFE         | 0.0217    | WSINDy/RFE       | 0/0 ✓     | 0.0048      | 0.0217     | 0.0222       |
| 5%    | pointwise     | spectral   | —        | 19/0 ✗            | 0.2287     | —                      | 0.5190    | —                | 2/0 ✗     | 0.2669      | 0.5437     | 2.3614       |
| 5%    | pointwise     | polynomial | 🔴        | 3/0 ✗             | 0.1123     | —                      | 0.5000    | PolyJet/STRidge  | 0/0 ✓     | 0.1254      | 0.5149     | 0.5325       |
| 5%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0063     | multigraph/MSTLS       | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0063      | 0.0536     | 0.0540       |
| 5%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0063     | mod_multigraph/MSTLS   | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0063      | 0.0536     | 0.0540       |
| 5%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0063     | multigraph/MSTLS       | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0063      | 0.0536     | 0.0540       |
| 5%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0063     | mod_multigraph/MSTLS   | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0063      | 0.0536     | 0.0540       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                   | fp/fn | E_coeff | Time (s) | vs. best  |
|:-------------------------------|:-----------------------------|:------|--------:|---------:|----------:|
| Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | —         |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0048  | 0.4      | +0.4%     |
| WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.3      | +0.4%     |
| PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.5      | +0.4%     |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0048  | 0.1      | +0.8%     |
| Multigraphs + MSTLS            | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.8%     |
| Modified Multigraphs + RFE     | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | +0.8%     |
| Multigraphs + RFE              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.8%     |
| WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +1.1%     |
| PolyJet + MSTLS                | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.8      | +1.1%     |
| PDE-FIND + MSTLS               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.3      | +1.1%     |
| PolyJet + RFE                  | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +1.1%     |
| PDE-FIND + RFE                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +1.1%     |
| WSINDy + MSTLS                 | polynomial, weak (IBP), m=12 | 0/0 ✓ | 2.8796  | 0.2      | +60074.4% |

### 2% noise

| Method                         | Best conf.                   | fp/fn | E_coeff | Time (s) | vs. best  |
|:-------------------------------|:-----------------------------|:------|--------:|---------:|----------:|
| Multigraphs + RFE              | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | —         |
| Modified Multigraphs + RFE     | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.0%     |
| Multigraphs + MSTLS            | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| Modified Multigraphs + MSTLS   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| WSINDy + RFE                   | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PDE-FIND + RFE                 | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PolyJet + RFE                  | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PDE-FIND + MSTLS               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.3      | +0.2%     |
| PolyJet + MSTLS                | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0048  | 0.7      | +0.2%     |
| Multigraphs + STRidge          | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.1      | +2.0%     |
| Modified Multigraphs + STRidge | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.1      | +2.0%     |
| PolyJet + STRidge              | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.4      | +2.4%     |
| PDE-FIND + STRidge             | polynomial, weak (IBP), m=12 | 0/0 ✓ | 0.0049  | 0.2      | +2.4%     |
| WSINDy + STRidge               | spectral, weak (IBP), m=12   | 0/0 ✓ | 0.0049  | 0.2      | +2.4%     |
| WSINDy + MSTLS                 | spectral, weak (IBP), m=12   | 0/0 ✓ | 2.8779  | 0.2      | +59740.9% |

### 5% noise

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|----------:|
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 1.1      | —         |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 0.3      | +0.0%     |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 0.2      | +0.3%     |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 0.1      | +0.3%     |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0065  | 0.0      | +3.5%     |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0065  | 0.1      | +3.5%     |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.2      | +3.9%     |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.5      | +3.9%     |
| PDE-FIND + STRidge             | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.2      | +3.9%     |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.1      | +4.2%     |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.2      | +4.2%     |
| WSINDy + MSTLS                 | spectral, weak (IBP), m=12      | 0/0 ✓ | 2.8737  | 0.2      | +45487.1% |
| WSINDy + RFE                   | spectral, pointwise             | 0/2 ✗ | 0.3407  | 0.0      | ✗         |
| PDE-FIND + RFE                 | spectral, pointwise             | 0/2 ✗ | 0.3407  | 0.1      | ✗         |
| PolyJet + RFE                  | spectral, pointwise             | 0/2 ✗ | 0.3407  | 0.1      | ✗         |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0084  | 0.0664 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 2/0   | 6        | 0.1340  | 0.3047 | 0.4      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 14/0  | 18       | 0.2579  | 0.5349 | 0.8      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0084  | 0.0664 | 0.4      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 2/0   | 6        | 0.1340  | 0.3047 | 0.4      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 16/0  | 12       | 0.2337  | 0.5249 | 0.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0084  | 0.0664 | 0.8      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 2/0   | 6        | 0.1340  | 0.3047 | 0.8      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 2/0   | 6        | 0.2669  | 0.5437 | 1.0      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 0/0   | 3        | 0.0084  | 0.0664 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 40/0  | 20       | 0.1497  | 0.2885 | 0.2      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 40/0  | 20       | 0.2197  | 0.5132 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 19/0  | 6        | 0.0501  | 0.0455 | 0.1      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 40/0  | 20       | 0.1497  | 0.2885 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 19/0  | 8        | 0.2287  | 0.5190 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0117  | 0.1133 | 0.8      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 4        | 0.0115  | 0.1133 | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 0/0   | 4        | 0.0278  | 0.2545 | 0.4      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 2/0   | 6        | 0.1007  | 0.4998 | 0.4      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0117  | 0.1133 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 4        | 0.0115  | 0.1133 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 0/0   | 4        | 0.0278  | 0.2545 | 0.3      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 0/0   | 4        | 0.1254  | 0.5149 | 0.3      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0117  | 0.1133 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 4        | 0.0115  | 0.1133 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 0/0   | 4        | 0.0278  | 0.2545 | 1.5      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 0/0   | 4        | 0.1254  | 0.5149 | 1.2      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 0/0   | 3        | 0.0116  | 0.1133 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0114  | 0.1133 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 0/0   | 3        | 0.0269  | 0.2545 | 0.2      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 3/0   | 4        | 0.1123  | 0.5000 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 0/0   | 3        | 0.0116  | 0.1133 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0114  | 0.1133 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 0/0   | 3        | 0.0269  | 0.2545 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 4/0   | 4        | 0.1327  | 0.5077 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.1      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 0.4      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 0.5      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0066  | 0.0536 | 0.4      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0063  | 0.0536 | 0.5      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3555  | 0.2544 | 0.1      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8796  | 2.8813 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.8796  | 2.8813 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8779  | 2.8793 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.8779  | 2.8793 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8737  | 2.8685 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.8737  | 2.8685 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3555  | 0.2544 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.8      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.8      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.4      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.7      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0049  | 0.0217 | 0.4      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.7      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0073  | 0.0770 | 0.6      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 1.1      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | yes | 0/0   | 4        | 0.0073  | 0.0770 | 0.6      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0063  | 0.0536 | 1.0      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3555  | 0.2544 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.4      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.5      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.4      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.5      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0065  | 0.0536 | 0.2      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 0.3      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0063  | 0.0536 | 0.3      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3555  | 0.2544 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8797  | 2.8812 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.8796  | 2.8813 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8779  | 2.8792 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.8779  | 2.8793 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0065  | 0.0536 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8737  | 2.8685 | 0.3      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.8737  | 2.8685 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3555  | 0.2544 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.6      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet          | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.6      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 1.5      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0049  | 0.0217 | 0.6      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 1.4      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0065  | 0.0536 | 0.5      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 1.1      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0066  | 0.0536 | 0.5      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0063  | 0.0536 | 1.1      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3555  | 0.2544 | 0.0      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0065  | 0.0770 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0065  | 0.0770 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0063  | 0.0536 | 0.2      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | yes | 0/0   | 3        | 0.0066  | 0.0770 | 0.1      | 🟢      |


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
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0061       | 1477.8  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0057       | 2682.2  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0061       | 1419.7  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0057       | 2504.0  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.0713       | 2330.4  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.0315       | 2821.4  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.0222       | 4037.5  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.0655       | 1731.0  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.0315       | 1425.0  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.0222       | 2331.3  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.0315       | 2841.3  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.0222       | 4037.9  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.0315       | 1413.3  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.0222       | 2329.2  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.1759       | 10382.8 |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.0773       | 13252.1 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.0540       | 19563.0 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.1578       | 2043.8  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.0773       | 1387.5  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.0540       | 1856.6  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.0773       | 13321.0 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.0540       | 19566.7 |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.0773       | 1381.1  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.0540       | 1856.5  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ac`.

| Source                                  | Wall time | CPU | Threads |
|:----------------------------------------|----------:|:----|:-------:|
| `sweep_paper_v2_ac_20260804_223951.csv` | ? min     | ?   | ?       |
