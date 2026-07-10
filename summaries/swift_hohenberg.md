# Swift-Hohenberg equation

The equation under consideration is

$$f_t = (r-1)f - 2\Delta f - \Delta^2 f - f^3, \qquad r = 0.3$$

which can be written in terms of $O(2)$-equivariant terms

| Operator     | True coefficient |
|:-------------|:-----------------|
| $f$          | $r-1 = -0.7$     |
| $\Delta f$   | $-2.0$           |
| $\Delta^2 f$ | $-1.0$           |
| $f^3$        | $-1.0$           |

or standard terms (PolyJet basis)

| Operator            | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------|:-----------|:---------|:-------|:---------|:-----------|
| $f$                 |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f^3$               |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xx} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |

Swift-Hohenberg contains a 4th-order biharmonic term $\Delta^2 f = f_{xxxx} + 2f_{xxyy} + f_{yyyy}$ which is the main source of difficulty.  All four libraries can represent every term individually, but the equivariant basis encodes the isotropic combination as a single operator, giving it a sparsity advantage

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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mg E_data | Std best         | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:----------|:-----------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | 🟩        | 0/0 ✓             | 0.0273     | 0.0209    | WSINDy/RFE       | 0/0 ✓     | 0.0274      | 0.0208     | 0.0197       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.2089     | 0.1511    | WSINDy/RFE       | 0/0 ✓     | 0.2110      | 0.1484     | 0.0914       |
| Clean | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0259     | 0.0138    | PDE-FIND/STRidge | 0/0 ✓     | 0.0257      | 0.0138     | 0.0108       |
| Clean | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0259     | 0.0139    | PDE-FIND/STRidge | 0/0 ✓     | 0.0258      | 0.0135     | 0.0109       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0224     | 0.0153    | PolyJet/STRidge  | 0/0 ✓     | 0.0230      | 0.0138     | 0.0110       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0224     | 0.0153    | PolyJet/STRidge  | 0/0 ✓     | 0.0230      | 0.0138     | 0.0110       |
| 2%    | pointwise     | spectral   | —        | 0/4 ✗             | 1.0000     | 0.9440    | —                | 0/5 ✗     | 1.0000      | 0.8727     | 0.6583       |
| 2%    | pointwise     | polynomial | —        | 0/4 ✗             | 1.0002     | 0.9436    | —                | 0/5 ✗     | 1.0002      | 0.9409     | 0.6010       |
| 2%    | weak (no IBP) | spectral   | 🟥        | 0/0 ✓             | 0.0096     | 0.1183    | PDE-FIND/MSTLS   | 0/0 ✓     | 0.0085      | 0.1183     | 0.1174       |
| 2%    | weak (no IBP) | polynomial | 🟥        | 0/0 ✓             | 0.0107     | 0.1180    | PolyJet/RFE      | 0/0 ✓     | 0.0097      | 0.1180     | 0.1166       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0079     | 0.1185    | PolyJet/RFE      | 0/0 ✓     | 0.0099      | 0.1184     | 0.1174       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0079     | 0.1185    | PolyJet/RFE      | 0/0 ✓     | 0.0099      | 0.1184     | 0.1169       |
| 5%    | pointwise     | spectral   | —        | 0/6 ✗             | 1.0000     | 0.9724    | —                | 0/5 ✗     | 1.0000      | 0.9404     | 0.8603       |
| 5%    | pointwise     | polynomial | —        | 0/4 ✗             | 1.0001     | 0.9723    | —                | 0/5 ✗     | 1.0000      | 0.9404     | 0.8442       |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0472     | 0.2848    | PolyJet/STRidge  | 0/0 ✓     | 0.0559      | 0.2845     | 0.2826       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0448     | 0.2839    | PolyJet/STRidge  | 0/0 ✓     | 0.0530      | 0.2837     | 0.2788       |
| 5%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0487     | 0.2847    | PDE-FIND/STRidge | 0/0 ✓     | 0.0548      | 0.2846     | 0.2825       |
| 5%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0487     | 0.2847    | PDE-FIND/STRidge | 0/0 ✓     | 0.0548      | 0.2846     | 0.2794       |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm   | d_tol/α   | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.005     | 0/0   | 7        | 0.0276  | 0.0208 | 0.2      | 🟢      |
| PDE-FIND    | Clean | RFE         |           | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.05      | 0/5   | 2        | 1.0000  | 0.8727 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 0/5   | 2        | 1.0000  | 0.8727 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.05      | 0/5   | 2        | 1.0000  | 0.9404 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 0/5   | 2        | 1.0000  | 0.9404 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.005     | 0/0   | 7        | 0.0276  | 0.0208 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 0/0   | 7        | 0.0274  | 0.0208 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.005     | 0/5   | 2        | 1.0000  | 0.9439 | 0.2      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 2/5   | 3        | 1.0001  | 0.7532 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.05      | 0/5   | 2        | 1.0000  | 0.9404 | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 6/5   | 4        | 1.0000  | 0.8806 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.05      | 0/0   | 7        | 0.0276  | 0.0208 | 0.5      | 🟢      |
| PolyJet     | Clean | RFE         |           | 0/0   | 7        | 0.0274  | 0.0208 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.005     | 0/5   | 2        | 1.0000  | 0.9439 | 0.3      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 1/5   | 3        | 1.0000  | 0.6951 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.01      | 0/6   | 1        | 1.0000  | 0.9724 | 0.3      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 1/5   | 3        | 1.0000  | 0.8872 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.1       | 0/0   | 4        | 0.0275  | 0.0209 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 0/0   | 4        | 0.0273  | 0.0209 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 0/0   | 4        | 0.0273  | 0.0209 | 18.6     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.005     | 0/4   | 2        | 1.0000  | 0.9440 | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 2/5   | 3        | 1.0000  | 0.6850 | 0.0      | 🔴      |
| Multigraphs | 2%    | Brute-force | L₀=1.0    | 0/6   | 1        | 1.0000  | 0.9441 | 17.9     | 🔴      |
| Multigraphs | 5%    | STRidge     | 0.005     | 0/6   | 1        | 1.0000  | 0.9724 | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 2/5   | 3        | 1.0000  | 0.8860 | 0.0      | 🔴      |
| Multigraphs | 5%    | Brute-force | L₀=1.0    | 0/6   | 1        | 1.0000  | 0.9724 | 17.9     | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm   | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.05    | 0/0   | 7        | 0.2122  | 0.1484 | 0.2      | 🟢      |
| PDE-FIND    | Clean | RFE         |         | 0/0   | 7        | 0.2110  | 0.1484 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.01    | 0/6   | 1        | 1.0001  | 0.9409 | 0.2      | 🔴      |
| PDE-FIND    | 2%    | RFE         |         | 1/6   | 2        | 1.0003  | 0.8649 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.1     | 0/7   | 0        | 1.0000  | 1.0000 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | RFE         |         | 1/6   | 2        | 1.0003  | 0.9391 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.05    | 0/0   | 7        | 0.2122  | 0.1484 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |         | 0/0   | 7        | 0.2110  | 0.1484 | 0.1      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.2     | 0/7   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| WSINDy      | 2%    | RFE         |         | 2/5   | 3        | 1.0002  | 0.7105 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.05    | 0/5   | 2        | 1.0000  | 0.9404 | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE         |         | 2/5   | 3        | 1.0001  | 0.8799 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.1     | 0/0   | 7        | 0.2122  | 0.1484 | 0.5      | 🟢      |
| PolyJet     | Clean | RFE         |         | 0/0   | 7        | 0.2110  | 0.1484 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.01    | 0/5   | 2        | 1.0002  | 0.9409 | 0.3      | 🔴      |
| PolyJet     | 2%    | RFE         |         | 2/5   | 4        | 1.0003  | 0.6451 | 0.1      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.2     | 0/7   | 0        | 1.0000  | 1.0000 | 0.3      | 🔴      |
| PolyJet     | 5%    | RFE         |         | 2/5   | 4        | 1.0002  | 0.8591 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.2     | 0/0   | 4        | 0.2089  | 0.1511 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |         | 2/1   | 4        | 0.2472  | 0.2349 | 0.0      | 🔴      |
| Multigraphs | Clean | Brute-force | L₀=1.0  | 0/5   | 1        | 0.9338  | 0.9302 | 18.7     | 🔴      |
| Multigraphs | 2%    | STRidge     | 0.005   | 0/4   | 2        | 1.0002  | 0.9436 | 0.1      | 🔴      |
| Multigraphs | 2%    | RFE         |         | 13/6  | 7        | 0.9997  | 0.6651 | 0.0      | 🔴      |
| Multigraphs | 2%    | Brute-force | L₀=1.0  | 0/5   | 1        | 1.0002  | 0.9436 | 17.6     | 🔴      |
| Multigraphs | 5%    | STRidge     | 0.005   | 0/4   | 2        | 1.0001  | 0.9723 | 0.1      | 🔴      |
| Multigraphs | 5%    | RFE         |         | 2/5   | 3        | 1.0000  | 0.8640 | 0.0      | 🔴      |
| Multigraphs | 5%    | Brute-force | L₀=1.0  | 0/6   | 1        | 1.0000  | 0.9724 | 17.6     | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm   | d_tol/α   | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 1.5       | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | no  | 2/0   | 9        | 0.0301  | 0.0158  | 1.5      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND    | Clean | STRidge     | 0.5       | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | yes | 0/7   | 0        | 1.0000  | 1.0000  | 1.6      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 1.0       | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.6      | 🟢      |
| PDE-FIND    | 2%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 1.0       | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.6      | 🟢      |
| PDE-FIND    | 2%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge     | 0.5       | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 12 | no  | 2/0   | 10       | 0.1000  | 0.2825  | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge     | 1.0       | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 12 | yes | 2/0   | 10       | 0.0992  | 0.2826  | 0.6      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.1       | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 8  | no  | 0/0   | 7        | 76.5564 | 70.8534 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.2      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 8  | yes | 0/0   | 7        | 74.2567 | 64.4902 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | no  | 6/0   | 9        | 6.0366  | 63.1279 | 0.4      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.5       | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | yes | 6/0   | 9        | 6.0511  | 63.2057 | 0.4      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 1.0       | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 8  | no  | 7/2   | 8        | 1.7297  | 41.0796 | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 1.0       | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 8  | yes | 7/2   | 8        | 1.7190  | 41.4895 | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge     | 4.0       | 12 | no  | 0/0   | 7        | 0.0257  | 0.0138  | 0.2      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 4  | no  | 5/0   | 12       | 0.0302  | 0.0156  | 3.8      | 🔴      |
| PolyJet     | Clean | RFE         |           | 12 | no  | 0/0   | 7        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge     | 3.0       | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.2      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 4  | yes | 0/7   | 0        | 1.0000  | 1.0000  | 2.8      | 🔴      |
| PolyJet     | Clean | RFE         |           | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 7        | 0.0095  | 0.1183  | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 4  | no  | 5/5   | 7        | 1.0002  | 0.4917  | 2.7      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0085  | 0.1183  | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 1.0       | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 4  | yes | 5/5   | 7        | 1.0002  | 0.4917  | 2.8      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.0      | 🟢      |
| PolyJet     | 5%    | STRidge     | 4.0       | 12 | no  | 0/0   | 7        | 0.0559  | 0.2845  | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 4  | no  | 2/5   | 4        | 1.0000  | 0.6174  | 2.8      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0584  | 0.2845  | 0.1      | 🟢      |
| PolyJet     | 5%    | STRidge     | 4.0       | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 4  | yes | 2/5   | 4        | 1.0000  | 0.6174  | 2.5      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge     | 1.5       | 12 | no  | 0/0   | 4        | 0.0259  | 0.0138  | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 4  | no  | 0/7   | 0        | 1.0000  | 1.0000  | 0.3      | 🔴      |
| Multigraphs | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 4        | 0.0260  | 0.0138  | 17.8     | 🟢      |
| Multigraphs | Clean | STRidge     | 1.5       | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | yes | 3/6   | 2        | 1.9990  | 0.8143  | 0.2      | 🔴      |
| Multigraphs | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 18.6     | 🟢      |
| Multigraphs | 2%    | STRidge     | 1.0       | 12 | no  | 0/0   | 4        | 0.0096  | 0.1183  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | no  | 0/7   | 0        | 1.0000  | 1.0000  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0100  | 0.1183  | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 4        | 0.0100  | 0.1183  | 18.0     | 🟢      |
| Multigraphs | 2%    | STRidge     | 1.0       | 12 | yes | 0/0   | 4        | 0.0079  | 0.1185  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | yes | 2/5   | 3        | 1.0002  | 0.5916  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0083  | 0.1185  | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 4        | 0.0083  | 0.1185  | 18.6     | 🟢      |
| Multigraphs | 5%    | STRidge     | 7.0       | 12 | no  | 0/0   | 4        | 0.0480  | 0.2848  | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | no  | 0/7   | 0        | 1.0000  | 1.0000  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0472  | 0.2848  | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 12 | no  | 0/0   | 4        | 0.0472  | 0.2848  | 18.7     | 🟢      |
| Multigraphs | 5%    | STRidge     | 4.0       | 12 | yes | 0/0   | 4        | 0.0496  | 0.2848  | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | yes | 2/5   | 3        | 1.0000  | 0.6086  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0487  | 0.2847  | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 12 | yes | 0/0   | 4        | 0.0487  | 0.2847  | 18.8     | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm   | d_tol/α   | m  | IBP | fp/fn | n_active | E_coeff | E_data  | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:---|:----|:------|:---------|:--------|:--------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.1       | 12 | no  | 0/0   | 9        | 0.0258  | 0.0135  | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | no  | 0/7   | 0        | 1.0000  | 1.0000  | 1.7      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 12 | no  | 0/0   | 7        | 0.0260  | 0.0138  | 0.0      | 🟢      |
| PDE-FIND    | Clean | STRidge     | 1.0       | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | yes | 0/7   | 0        | 1.0000  | 1.0000  | 1.7      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 1.5       | 12 | no  | 0/0   | 7        | 0.0108  | 0.1180  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 8  | no  | 1/4   | 6        | 0.9468  | 0.4655  | 1.0      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0097  | 0.1180  | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.5       | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 8  | yes | 1/4   | 6        | 0.9469  | 0.4655  | 1.0      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge     | 1.5       | 12 | no  | 0/0   | 7        | 0.0530  | 0.2837  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 8  | no  | 1/5   | 3        | 1.0518  | 0.6273  | 1.1      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0556  | 0.2837  | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge     | 1.5       | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 8  | yes | 1/5   | 3        | 1.0511  | 0.6276  | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.2       | 12 | no  | 0/0   | 7        | 0.0258  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 8  | no  | 0/0   | 7        | 76.6159 | 70.8725 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | no  | 0/0   | 7        | 0.0260  | 0.0138  | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 8  | yes | 0/0   | 7        | 74.2567 | 64.4902 | 0.7      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 7        | 0.0108  | 0.1180  | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | no  | 6/0   | 9        | 6.0434  | 63.1179 | 0.4      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0097  | 0.1180  | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.5       | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | yes | 6/0   | 9        | 6.0508  | 63.2065 | 0.5      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 1.0       | 12 | no  | 0/0   | 7        | 0.0530  | 0.2837  | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 8  | no  | 11/1  | 10       | 1.0841  | 56.3276 | 0.7      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0556  | 0.2837  | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 1.0       | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 8  | yes | 9/2   | 8        | 1.8304  | 44.3223 | 0.6      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge     | 3.0       | 12 | no  | 0/0   | 7        | 0.0258  | 0.0138  | 0.2      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 8  | no  | 4/0   | 11       | 0.1316  | 0.0147  | 2.2      | 🔴      |
| PolyJet     | Clean | RFE         |           | 12 | no  | 0/0   | 7        | 0.0260  | 0.0138  | 0.1      | 🟢      |
| PolyJet     | Clean | STRidge     | 3.0       | 12 | yes | 0/0   | 7        | 0.0230  | 0.0138  | 0.2      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 4  | yes | 0/7   | 0        | 1.0000  | 1.0000  | 3.1      | 🔴      |
| PolyJet     | Clean | RFE         |           | 12 | yes | 0/0   | 7        | 0.0233  | 0.0138  | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge     | 2.0       | 12 | no  | 0/0   | 7        | 0.0108  | 0.1180  | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 4  | no  | 14/5  | 16       | 1.0003  | 0.4052  | 2.9      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0097  | 0.1180  | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 1.5       | 12 | yes | 0/0   | 7        | 0.0109  | 0.1184  | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 8  | yes | 10/3  | 14       | 0.9283  | 0.4140  | 2.1      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0099  | 0.1184  | 0.0      | 🟢      |
| PolyJet     | 5%    | STRidge     | 2.0       | 12 | no  | 0/0   | 7        | 0.0530  | 0.2837  | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 4  | no  | 9/6   | 10       | 1.0002  | 0.5579  | 3.0      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 12 | no  | 0/0   | 7        | 0.0556  | 0.2837  | 0.0      | 🟢      |
| PolyJet     | 5%    | STRidge     | 3.0       | 12 | yes | 0/0   | 7        | 0.0548  | 0.2846  | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 4  | yes | 9/6   | 10       | 1.0002  | 0.5518  | 3.2      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 12 | yes | 0/0   | 7        | 0.0573  | 0.2846  | 0.0      | 🟢      |
| Multigraphs | Clean | STRidge     | 3.0       | 12 | no  | 0/0   | 4        | 0.0259  | 0.0139  | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 8  | no  | 0/7   | 0        | 1.0000  | 1.0000  | 0.2      | 🔴      |
| Multigraphs | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0261  | 0.0139  | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 4        | 0.0261  | 0.0139  | 20.3     | 🟢      |
| Multigraphs | Clean | STRidge     | 1.5       | 12 | yes | 0/0   | 4        | 0.0224  | 0.0153  | 0.1      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | yes | 3/6   | 2        | 1.9974  | 0.8147  | 0.2      | 🔴      |
| Multigraphs | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 4        | 0.0225  | 0.0153  | 18.3     | 🟢      |
| Multigraphs | 2%    | STRidge     | 2.0       | 12 | no  | 0/0   | 4        | 0.0107  | 0.1180  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 8  | no  | 0/7   | 0        | 1.0000  | 1.0000  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0110  | 0.1180  | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 4        | 0.0110  | 0.1180  | 17.9     | 🟢      |
| Multigraphs | 2%    | STRidge     | 4.0       | 12 | yes | 0/0   | 4        | 0.0079  | 0.1185  | 0.1      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | yes | 8/6   | 3        | 1.0007  | 0.6560  | 0.2      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0083  | 0.1185  | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.01   | 12 | yes | 0/0   | 4        | 0.0083  | 0.1185  | 18.5     | 🟢      |
| Multigraphs | 5%    | STRidge     | 2.0       | 12 | no  | 0/0   | 4        | 0.0455  | 0.2840  | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | no  | 7/5   | 4        | 1.0006  | 0.5725  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0448  | 0.2839  | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 12 | no  | 0/0   | 4        | 0.0448  | 0.2839  | 17.8     | 🟢      |
| Multigraphs | 5%    | STRidge     | 2.0       | 12 | yes | 0/0   | 4        | 0.0496  | 0.2848  | 0.1      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | yes | 7/5   | 4        | 1.0006  | 0.5725  | 0.2      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0487  | 0.2847  | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 12 | yes | 0/0   | 4        | 0.0487  | 0.2847  | 18.7     | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0149       | 10372.4 |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0108       | 15893.5 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0109       | 28641.5 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0538       | 10430.3 |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0123       | 15952.5 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0109       | 28643.3 |
| Clean | weak (IBP)    | spectral   | 4  | yes | 0.2767       | 10380.5 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0359       | 15897.9 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0110       | 28695.1 |
| Clean | weak (IBP)    | polynomial | 4  | yes | 0.2761       | 10448.2 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0358       | 15962.4 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0110       | 28691.7 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.4727       | 692.9   |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.2980       | 1432.1  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.1174       | 3350.2  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.4397       | 8016.7  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.2704       | 13410.6 |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.1166       | 25224.6 |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 0.4723       | 691.0   |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.2972       | 1432.2  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.1174       | 3350.1  |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 0.4446       | 8011.6  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.2907       | 13414.0 |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.1169       | 25249.1 |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.5978       | 744.2   |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.4616       | 1578.6  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.2826       | 4258.3  |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.5392       | 3051.7  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.4403       | 6492.0  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.2788       | 13259.2 |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 0.5977       | 741.0   |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.4607       | 1579.0  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.2825       | 4258.1  |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 0.5420       | 3037.9  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.4496       | 6488.1  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.2794       | 13263.5 |
---
| Source                                     | Wall time | CPU                                      | Threads |
|:-------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260710_002834_sh_mstls_fixed.csv` | 197.5 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 4       |
