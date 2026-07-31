# Cahn-Hilliard equation

The equation under consideration is

$$f_t = -\Delta^2 f + \Delta(f^3 - f) = -\Delta f + 3f^2\Delta f + 6f\|\nabla f\|^2 - \varepsilon^2\Delta^2 f, \qquad \varepsilon = 0.2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator               | True coefficient         |
|:-----------------------|:-------------------------|
| $\Delta f$             | $-1.0$                   |
| $f^2\cdot\Delta f$     | $+3.0$                   |
| $f\cdot\|\nabla f\|^2$ | $+6.0$                   |
| $\Delta^2 f$           | $-\varepsilon^2 = -0.04$ |

or standard terms (PolyJet basis)

| Operator                  | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------------|:-----------|:---------|:-------|:---------|:-----------|
| $\partial_{xx} f$         |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$         |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^2\cdot\partial_{xx} f$ |            | 🟢        | 🔴      | 🟢        | 🟢          |
| $f^2\cdot\partial_{yy} f$ |            | 🟢        | 🔴      | 🟢        | 🟢          |
| $f\cdot f_x^2$            |            | 🔴        | 🔴      | 🟢        | 🟢          |
| $f\cdot f_y^2$            |            | 🔴        | 🔴      | 🟢        | 🟢          |
| $\partial_{xxxx} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$       |            | 🟢        | 🟢      | 🟢        | 🟢          |

Cahn-Hilliard combines second and fourth-order terms with nonlinear products.  PDE-FIND cannot represent $f\cdot f_x^2$ (products of derivatives).  WSINDy cannot represent either $f^2\cdot f_{xx}$ or $f\cdot f_x^2$ individually, but the combination $3f^2\Delta f + 6f\|\nabla f\|^2$ equals $\Delta(f^3)$, which IS representable as $(f^3)_{xx} + (f^3)_{yy}$ in the WSINDy basis.

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
| Clean | pointwise     | spectral   | 🟥        | 0/0 ✓             | 0.1221     | multigraph/STRidge | 0.2350    | WSINDy/STRidge | 0/0 ✓     | 0.1200      | 0.2360     | 0.2711       |
| Clean | pointwise     | polynomial | —        | 35/0 ✗            | 0.9960     | —                  | 0.6740    | —              | 10/0 ✗    | 0.6951      | 0.5600     | 13.7157      |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0419     | multigraph/RFE     | 0.0578    | WSINDy/RFE     | 0/0 ✓     | 0.0400      | 0.0597     | 0.0713       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0553     | multigraph/RFE     | 0.0853    | WSINDy/RFE     | 0/0 ✓     | 0.0389      | 0.0604     | 0.1714       |
| Clean | weak (IBP)    | spectral   | 🟥        | 0/0 ✓             | 0.0522     | multigraph/RFE     | 0.0789    | WSINDy/RFE     | 0/0 ✓     | 0.0512      | 0.0798     | 0.0958       |
| Clean | weak (IBP)    | polynomial | 🟥        | 0/0 ✓             | 0.0717     | multigraph/RFE     | 0.1280    | PolyJet/RFE    | 0/0 ✓     | 0.0664      | 0.1132     | 0.2349       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                | Best conf.                      | fp/fn  | E_coeff  | Time (s) | vs. best |
|:----------------------|:--------------------------------|:-------|---------:|---------:|---------:|
| WSINDy + RFE          | polynomial, weak (no IBP), m=12 | 0/0 ✓  | 0.0389   | 0.6      | —        |
| Multigraphs + RFE     | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0419   | 0.4      | +7.8%    |
| Multigraphs + MSTLS   | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0419   | 2.3      | +7.8%    |
| Multigraphs + STRidge | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0424   | 1.9      | +9.1%    |
| PolyJet + RFE         | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0434   | 1.0      | +11.7%   |
| PolyJet + MSTLS       | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0434   | 37.8     | +11.7%   |
| PolyJet + STRidge     | spectral, weak (no IBP), m=12   | 0/0 ✓  | 0.0442   | 11.8     | +13.7%   |
| WSINDy + STRidge      | spectral, weak (no IBP), m=4    | 0/0 ✓  | 0.0832   | 4.3      | +114.1%  |
| WSINDy + MSTLS        | polynomial, weak (no IBP), m=12 | 2/0 ✗  | 239.0799 | 4.6      | ✗        |
| PDE-FIND + STRidge    | spectral, weak (no IBP), m=12   | 5/2 ✗  | 4.2677   | 6.4      | ✗        |
| PDE-FIND + MSTLS      | polynomial, weak (no IBP), m=8  | 25/2 ✗ | 1.0162   | 8.2      | ✗        |
| PDE-FIND + RFE        | spectral, pointwise             | 2/7 ✗  | 0.8004   | 0.7      | ✗        |

### 2% noise

| Method | Best conf. | fp/fn | E_coeff | Time (s) | vs. best |
|:-------|:-----------|:------|--------:|---------:|---------:|

### 5% noise

| Method | Best conf. | fp/fn | E_coeff | Time (s) | vs. best |
|:-------|:-----------|:------|--------:|---------:|---------:|

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 30/2  | 48       | 0.8855  | 0.7491 | 7.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/7   | 4        | 0.8004  | 0.7664 | 0.7      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 0/0   | 7        | 0.1200  | 0.2360 | 4.1      | 🟢      |
| WSINDy      | Clean | RFE       |         | 0/0   | 7        | 0.1205  | 0.2360 | 0.6      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 0/0   | 9        | 0.1226  | 0.2348 | 12.7     | 🟢      |
| PolyJet     | Clean | RFE       |         | 0/0   | 9        | 0.1235  | 0.2348 | 1.0      | 🟢      |
| Multigraphs | Clean | STRidge   | 1.0     | 0/0   | 4        | 0.1221  | 0.2350 | 2.4      | 🟢      |
| Multigraphs | Clean | RFE       |         | 0/0   | 4        | 0.1231  | 0.2350 | 0.5      | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 0.01    | 37/2  | 56       | 0.9984  | 0.7513 | 7.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 2/7   | 4        | 0.9892  | 0.7796 | 0.7      | 🔴      |
| WSINDy      | Clean | STRidge   | 1.0     | 10/0  | 10       | 0.6956  | 0.5600 | 4.4      | 🔴      |
| WSINDy      | Clean | RFE       |         | 10/0  | 10       | 0.6951  | 0.5600 | 0.6      | 🔴      |
| PolyJet     | Clean | STRidge   | 0.1     | 23/0  | 32       | 0.9962  | 0.6084 | 15.1     | 🔴      |
| PolyJet     | Clean | RFE       |         | 3/7   | 5        | 0.9917  | 0.8202 | 1.0      | 🔴      |
| Multigraphs | Clean | STRidge   | 0.001   | 35/0  | 20       | 0.9960  | 0.6740 | 2.2      | 🔴      |
| Multigraphs | Clean | RFE       |         | 6/7   | 4        | 1.0023  | 0.7429 | 0.5      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff   | E_data    | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:----------|:----------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 10.0    | 12 | no  | 5/2   | 12       | 4.2677    | 0.6828    | 6.4      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | no  | 32/2  | 49       | 4.1336    | 0.6586    | 8.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 4  | no  | 2/7   | 4        | 0.8013    | 0.7504    | 0.8      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 1.0     | 8  | yes | 10/2  | 20       | 0.8636    | 0.7227    | 5.8      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 30/2  | 46       | 1.4703    | 0.6853    | 7.9      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/7   | 4        | 0.8010    | 0.7347    | 0.7      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 4  | no  | 0/0   | 7        | 0.0832    | 0.1575    | 4.3      | 🟢      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 2/0   | 9        | 3393.2545 | 609.0389  | 4.6      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0400    | 0.0597    | 0.8      | 🟢      |
| WSINDy      | Clean | STRidge   | 1.0     | 12 | yes | 2/0   | 9        | 0.0637    | 0.0753    | 4.7      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 9        | 6093.0620 | 1975.2800 | 4.9      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0512    | 0.0798    | 0.7      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 9        | 0.0442    | 0.0578    | 11.8     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 9        | 0.0434    | 0.0578    | 37.8     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 9        | 0.0434    | 0.0578    | 1.0      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 9        | 0.0535    | 0.0789    | 12.9     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 9        | 0.0525    | 0.0789    | 40.7     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 9        | 0.0525    | 0.0789    | 1.1      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 4        | 0.0424    | 0.0578    | 1.9      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 4        | 0.0419    | 0.0578    | 2.3      | 🟢      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0419    | 0.0578    | 0.4      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 0/0   | 4        | 0.0532    | 0.0789    | 2.4      | 🟢      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 6        | 0.0625    | 0.0751    | 2.3      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0522    | 0.0789    | 0.5      | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff  | E_data    | Time (s) | Exact? |
|:------------|:------|:----------|:--------|:---|:----|:------|:---------|:---------|:----------|:---------|:-------|
| PDE-FIND    | Clean | STRidge   | 10.0    | 12 | no  | 8/2   | 18       | 1.1857   | 0.6756    | 5.9      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 8  | no  | 25/2  | 40       | 1.0162   | 0.7034    | 8.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | no  | 2/7   | 4        | 0.9885   | 0.7297    | 0.7      | 🔴      |
| PDE-FIND    | Clean | STRidge   | 10.0    | 12 | yes | 7/2   | 16       | 1.0163   | 0.7005    | 5.8      | 🔴      |
| PDE-FIND    | Clean | MSTLS     | 0.75    | 12 | yes | 30/2  | 46       | 1.0072   | 0.6853    | 8.2      | 🔴      |
| PDE-FIND    | Clean | RFE       |         | 8  | yes | 2/7   | 4        | 0.9880   | 0.7347    | 0.7      | 🔴      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | no  | 2/0   | 9        | 0.0588   | 0.0570    | 4.4      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | no  | 2/0   | 9        | 239.0799 | 607.9449  | 4.6      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0389   | 0.0604    | 0.6      | 🟢      |
| WSINDy      | Clean | STRidge   | 10.0    | 12 | yes | 2/0   | 9        | 0.1046   | 0.1630    | 4.4      | 🔴      |
| WSINDy      | Clean | MSTLS     | 0.25    | 12 | yes | 23/0  | 16       | 413.1106 | 1926.2061 | 5.2      | 🔴      |
| WSINDy      | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0672   | 0.1713    | 0.6      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 9        | 0.0527   | 0.0744    | 12.3     | 🟢      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 9        | 0.0521   | 0.0744    | 40.8     | 🟢      |
| PolyJet     | Clean | RFE       |         | 12 | no  | 0/0   | 9        | 0.0521   | 0.0744    | 1.0      | 🟢      |
| PolyJet     | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 10       | 0.0926   | 0.0981    | 12.7     | 🔴      |
| PolyJet     | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 11       | 0.0871   | 0.1096    | 40.0     | 🔴      |
| PolyJet     | Clean | RFE       |         | 12 | yes | 0/0   | 9        | 0.0664   | 0.1132    | 1.2      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 5/0   | 7        | 0.1418   | 0.0704    | 2.2      | 🔴      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | no  | 2/0   | 6        | 0.0637   | 0.0846    | 2.5      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 4        | 0.0553   | 0.0853    | 0.5      | 🟢      |
| Multigraphs | Clean | STRidge   | 10.0    | 8  | yes | 2/0   | 6        | 0.2591   | 0.3084    | 2.1      | 🔴      |
| Multigraphs | Clean | MSTLS     | 0.25    | 12 | yes | 2/0   | 6        | 0.0932   | 0.1244    | 2.4      | 🔴      |
| Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 4        | 0.0717   | 0.1280    | 0.5      | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.1807       | 11574.8 |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.1025       | 7441.2  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0713       | 7501.5  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 4.7437       | 8970.8  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.4406       | 6618.7  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.1714       | 6893.7  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.1490       | 9562.9  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0958       | 7570.0  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.4561       | 8737.4  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.2349       | 6874.1  |
---
## Provenance

Generated by `scripts/build_summary.py` from the sweep CSV(s) below. Regenerate
after code or sweep changes with `python scripts/build_summary.py <csv> --equation ch`.

| Source                              | Wall time | CPU | Threads |
|:------------------------------------|----------:|:----|:-------:|
| `sweep_20260730_045148_ch_full.csv` | ? min     | ?   | ?       |
