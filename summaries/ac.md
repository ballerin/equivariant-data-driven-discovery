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
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0080     | multigraph/RFE     | 0.0664    | PolyJet/RFE      | 0/0 ✓     | 0.0080      | 0.0664     | —            |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0113     | multigraph/RFE     | 0.1133    | PolyJet/RFE      | 0/0 ✓     | 0.0114      | 0.1133     | —            |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge | 0.0032    | PolyJet/STRidge  | 0/0 ✓     | 0.0048      | 0.0032     | —            |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/STRidge | 0.0032    | PolyJet/STRidge  | 0/0 ✓     | 0.0048      | 0.0032     | —            |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0050     | multigraph/STRidge | 0.0034    | PDE-FIND/STRidge | 0/0 ✓     | 0.0050      | 0.0034     | —            |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0050     | multigraph/STRidge | 0.0034    | PDE-FIND/STRidge | 0/0 ✓     | 0.0050      | 0.0034     | —            |
| 2%    | pointwise     | spectral   | —        | 40/0 ✗            | 0.1497     | —                  | 0.2885    | —                | 2/0 ✗     | 0.1340      | 0.3047     | —            |
| 2%    | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0265     | multigraph/STRidge | 0.2545    | PDE-FIND/STRidge | 0/0 ✓     | 0.0273      | 0.2545     | —            |
| 2%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0048     | multigraph/RFE     | 0.0217    | WSINDy/RFE       | 0/0 ✓     | 0.0048      | 0.0217     | —            |
| 2%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048     | multigraph/RFE     | 0.0217    | WSINDy/RFE       | 0/0 ✓     | 0.0048      | 0.0217     | —            |
| 2%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0051     | multigraph/RFE     | 0.0287    | WSINDy/RFE       | 0/0 ✓     | 0.0051      | 0.0287     | —            |
| 2%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0051     | multigraph/RFE     | 0.0287    | WSINDy/RFE       | 0/0 ✓     | 0.0051      | 0.0287     | —            |
| 5%    | pointwise     | spectral   | —        | 5/0 ✗             | 0.2357     | —                  | 0.5288    | —                | 2/0 ✗     | 0.2669      | 0.5437     | —            |
| 5%    | pointwise     | polynomial | 🔴        | 4/0 ✗             | 0.1299     | —                  | 0.5077    | WSINDy/STRidge   | 0/0 ✓     | 0.1222      | 0.5149     | —            |
| 5%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0063     | multigraph/MSTLS   | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0063      | 0.0536     | —            |
| 5%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0063     | multigraph/MSTLS   | 0.0536    | PolyJet/MSTLS    | 0/0 ✓     | 0.0063      | 0.0536     | —            |
| 5%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0069     | multigraph/RFE     | 0.0713    | PDE-FIND/MSTLS   | 0/0 ✓     | 0.0069      | 0.0713     | —            |
| 5%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0069     | multigraph/RFE     | 0.0713    | PDE-FIND/MSTLS   | 0/0 ✓     | 0.0069      | 0.0713     | —            |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|----------:|
| Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | —         |
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.4      | +0.4%     |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.2      | +0.4%     |
| WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.2      | +0.4%     |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.8%     |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.8%     |
| Multigraphs + RFE              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.8%     |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.8%     |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.7      | +1.1%     |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.3      | +1.1%     |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +1.1%     |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +1.1%     |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +1.1%     |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 2.8299  | 0.2      | +59028.8% |

### 2% noise

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|----------:|
| Multigraphs + RFE              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | —         |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.0%     |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.0%     |
| WSINDy + RFE                   | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.1      | +0.2%     |
| PDE-FIND + RFE                 | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.0      | +0.2%     |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.8      | +0.2%     |
| PDE-FIND + MSTLS               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0048  | 0.3      | +0.2%     |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0049  | 0.1      | +2.0%     |
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0049  | 0.1      | +2.0%     |
| PolyJet + STRidge              | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0049  | 0.5      | +2.4%     |
| PDE-FIND + STRidge             | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0049  | 0.2      | +2.4%     |
| WSINDy + STRidge               | spectral, weak (no IBP), m=12   | 0/0 ✓ | 0.0049  | 0.2      | +2.4%     |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 2.8283  | 0.2      | +58707.4% |

### 5% noise

| Method                         | Best conf.                      | fp/fn | E_coeff | Time (s) | vs. best  |
|:-------------------------------|:--------------------------------|:------|--------:|---------:|----------:|
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 1.0      | —         |
| PDE-FIND + MSTLS               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 0.4      | +0.0%     |
| Multigraphs + MSTLS            | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 0.1      | +0.3%     |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0063  | 0.1      | +0.3%     |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0065  | 0.0      | +3.4%     |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=8  | 0/0 ✓ | 0.0065  | 0.0      | +3.4%     |
| PDE-FIND + STRidge             | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.2      | +4.0%     |
| WSINDy + STRidge               | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.2      | +4.0%     |
| PolyJet + STRidge              | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0065  | 0.4      | +4.0%     |
| Multigraphs + STRidge          | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.1      | +4.3%     |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 0.0066  | 0.1      | +4.3%     |
| WSINDy + MSTLS                 | polynomial, weak (no IBP), m=12 | 0/0 ✓ | 2.8245  | 0.2      | +44779.7% |
| PDE-FIND + RFE                 | polynomial, weak (no IBP), m=4  | 0/2 ✗ | 0.3299  | 0.0      | ✗         |
| WSINDy + RFE                   | polynomial, weak (no IBP), m=4  | 0/2 ✗ | 0.3299  | 0.0      | ✗         |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=4  | 0/2 ✗ | 0.3299  | 0.1      | ✗         |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0084  | 0.0664 | 0.4      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 2/0   | 6        | 0.1340  | 0.3047 | 0.5      | 🔴      |
| PDE-FIND         | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 14/0  | 18       | 0.2579  | 0.5349 | 0.6      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0084  | 0.0664 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 2/0   | 6        | 0.1340  | 0.3047 | 0.3      | 🔴      |
| WSINDy           | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 16/0  | 12       | 0.2337  | 0.5249 | 0.4      | 🔴      |
| WSINDy           | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0084  | 0.0664 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 2/0   | 6        | 0.1340  | 0.3047 | 0.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 2/0   | 6        | 0.2669  | 0.5437 | 1.0      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 0/0   | 3        | 0.0084  | 0.0664 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.001   | 40/0  | 20       | 0.1497  | 0.2885 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 5/0   | 5        | 0.2357  | 0.5288 | 0.2      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 11/0  | 6        | 0.0503  | 0.0453 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.001   | 40/0  | 20       | 0.1497  | 0.2885 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 0/2   | 2        | 0.3397  | 0.4072 | 0.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 19/0  | 8        | 0.2287  | 0.5190 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 0/2   | 2        | 0.3407  | 0.5669 | 0.0      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0116  | 0.1133 | 0.4      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 0/0   | 4        | 0.0114  | 0.1133 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.01    | 0/0   | 4        | 0.0273  | 0.2545 | 0.4      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 0/2   | 2        | 0.3291  | 0.4072 | 0.1      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.01    | 2/0   | 6        | 0.0792  | 0.4998 | 0.5      | 🔴      |
| PDE-FIND         | 5%    | RFE       |         | 0/2   | 2        | 0.3304  | 0.5669 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0116  | 0.1133 | 0.3      | 🟢      |
| WSINDy           | Clean | RFE       |         | 0/0   | 4        | 0.0114  | 0.1133 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 0/0   | 4        | 0.0273  | 0.2545 | 0.3      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 0/2   | 2        | 0.3291  | 0.4072 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.01    | 0/0   | 4        | 0.1222  | 0.5149 | 0.3      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 0/2   | 2        | 0.3304  | 0.5669 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.01    | 0/0   | 4        | 0.0116  | 0.1133 | 0.8      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 4        | 0.0114  | 0.1133 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.01    | 0/0   | 4        | 0.0273  | 0.2545 | 0.8      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 0/2   | 2        | 0.3291  | 0.4072 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.01    | 0/0   | 4        | 0.1222  | 0.5149 | 0.9      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 0/2   | 2        | 0.3304  | 0.5669 | 0.1      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.01    | 0/0   | 3        | 0.0115  | 0.1133 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0113  | 0.1133 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 0/0   | 3        | 0.0265  | 0.2545 | 0.2      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 0/2   | 2        | 0.3291  | 0.4072 | 0.0      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 40/0  | 20       | 0.0793  | 0.4923 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 0/2   | 2        | 0.3304  | 0.5669 | 0.0      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 0/0   | 3        | 0.0115  | 0.1133 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0113  | 0.1133 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 0/0   | 3        | 0.0265  | 0.2545 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 0/2   | 2        | 0.3291  | 0.4072 | 0.0      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.01    | 4/0   | 4        | 0.1299  | 0.5077 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 0/2   | 2        | 0.3304  | 0.5669 | 0.0      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0052  | 0.0287 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 0.3      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 0.4      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0074  | 0.0713 | 0.3      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0069  | 0.0713 | 0.4      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3450  | 0.3001 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8796  | 2.8813 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0766  | 2.7381 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8779  | 2.8793 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | yes | 0/0   | 4        | 0.0052  | 0.0287 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0735  | 2.7348 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0066  | 0.0536 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8737  | 2.8685 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0074  | 0.0713 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0628  | 2.7195 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3450  | 0.3001 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet          | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.5      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.8      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0052  | 0.0287 | 0.5      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.8      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 0/0   | 4        | 0.0073  | 0.0770 | 0.6      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 0.9      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3410  | 0.3345 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 8  | yes | 0/0   | 4        | 0.0155  | 0.1216 | 0.6      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0069  | 0.0713 | 0.9      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3450  | 0.3001 | 0.0      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs      | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0052  | 0.0287 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 10.0    | 12 | yes | 0/0   | 3        | 0.0074  | 0.0713 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0052  | 0.0287 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0066  | 0.0770 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0074  | 0.0713 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.0      | 🟢      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND         | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.2      | 🟢      |
| PDE-FIND         | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.3      | 🟢      |
| PDE-FIND         | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND         | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0052  | 0.0287 | 0.2      | 🟢      |
| PDE-FIND         | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.3      | 🟢      |
| PDE-FIND         | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0065  | 0.0536 | 0.2      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 0.4      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3299  | 0.3345 | 0.0      | 🔴      |
| PDE-FIND         | 5%    | STRidge   | 1.0     | 12 | yes | 0/0   | 4        | 0.0074  | 0.0713 | 0.2      | 🟢      |
| PDE-FIND         | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0069  | 0.0713 | 0.3      | 🟢      |
| PDE-FIND         | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3336  | 0.3001 | 0.0      | 🔴      |
| WSINDy           | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8299  | 2.8812 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.0      | 🟢      |
| WSINDy           | Clean | STRidge   | 0.01    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.2      | 🟢      |
| WSINDy           | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 3.0014  | 2.7381 | 0.2      | 🟢      |
| WSINDy           | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8283  | 2.8792 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| WSINDy           | 2%    | STRidge   | 0.01    | 12 | yes | 0/0   | 4        | 0.0052  | 0.0287 | 0.2      | 🟢      |
| WSINDy           | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.9985  | 2.7348 | 0.2      | 🟢      |
| WSINDy           | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| WSINDy           | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0065  | 0.0536 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 2.8245  | 2.8685 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3299  | 0.3345 | 0.0      | 🔴      |
| WSINDy           | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0074  | 0.0713 | 0.2      | 🟢      |
| WSINDy           | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 2.9883  | 2.7195 | 0.2      | 🟢      |
| WSINDy           | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3336  | 0.3001 | 0.0      | 🔴      |
| PolyJet          | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.5      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0049  | 0.0032 | 0.0      | 🟢      |
| PolyJet          | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.4      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.7      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| PolyJet          | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.4      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.8      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/0   | 4        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0052  | 0.0287 | 0.5      | 🟢      |
| PolyJet          | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.7      | 🟢      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/0   | 4        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| PolyJet          | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 4        | 0.0065  | 0.0536 | 0.4      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0063  | 0.0536 | 1.0      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 0/2   | 2        | 0.3299  | 0.3345 | 0.1      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 4        | 0.0074  | 0.0713 | 0.5      | 🟢      |
| PolyJet          | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 4        | 0.0069  | 0.0713 | 0.9      | 🟢      |
| PolyJet          | 5%    | RFE       |         | 8  | yes | 0/2   | 2        | 0.3336  | 0.3001 | 0.0      | 🔴      |
| Multigraphs      | Clean | STRidge   | 0.1     | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs      | Clean | STRidge   | 1.0     | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.01    | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0052  | 0.0287 | 0.1      | 🟢      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.1      | 🟢      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0065  | 0.0770 | 0.0      | 🟢      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0074  | 0.0713 | 0.1      | 🟢      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.1      | 🟢      |
| Multigraphs      | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.01    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.1      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0050  | 0.0034 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/0   | 3        | 0.0048  | 0.0217 | 0.0      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0052  | 0.0287 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.1      | 🟢      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0051  | 0.0287 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | no  | 0/0   | 3        | 0.0066  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 0/0   | 3        | 0.0063  | 0.0536 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 8  | no  | 0/0   | 3        | 0.0065  | 0.0770 | 0.0      | 🟢      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 12 | yes | 0/0   | 3        | 0.0074  | 0.0713 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.1      | 🟢      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/0   | 3        | 0.0069  | 0.0713 | 0.0      | 🟢      |

---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ac`.

| Source                                  | Wall time | CPU                                      | Threads |
|:----------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260801_223445_ac_full.csv`     | 33.8 min  | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
| `sweep_20260802_142233_mod_mg_full.csv` | 100.3 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
