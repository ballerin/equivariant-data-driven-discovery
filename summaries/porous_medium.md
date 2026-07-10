# Porous Medium equation

The equation under consideration is

$$f_t = \nabla\cdot(f^2\nabla f) = f^2\Delta f + 2f\|\nabla f\|^2$$

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

PDE-FIND cannot represent $f\cdot f_x^2$ or $f\cdot f_y^2$ (products of derivatives), so its minimum fn is 2.  WSINDy represents PME exactly via the chain rule: $f_t = \frac13\partial_x^2(f^3) + \frac13\partial_y^2(f^3)$, which expands to all 4 jet terms.

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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mg E_data | Std best | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:----------|:---------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0011     | 0.0019    | WSINDy   | 0/0 ✓     | 0.0010      | 0.0019     | 0.0009       |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0008     | 0.0051    | PolyJet  | 0/0 ✓     | 0.0008      | 0.0051     | 0.0024       |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0008     | 0.0006    | WSINDy   | 0/0 ✓     | 0.0007      | 0.0007     | 0.0000       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0007     | 0.0007    | WSINDy   | 0/0 ✓     | 0.0007      | 0.0007     | 0.0001       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0008     | 0.0006    | WSINDy   | 0/0 ✓     | 0.0007      | 0.0007     | 0.0000       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0007     | 0.0007    | WSINDy   | 0/0 ✓     | 0.0007      | 0.0007     | 0.0001       |
| 2%    | pointwise     | spectral   | —        | 0/4 ✗             | 1.0000     | 1.0000    | —        | 0/4 ✗     | 1.0000      | 1.0000     | 0.8522       |
| 2%    | pointwise     | polynomial | —        | 0/4 ✗             | 1.0000     | 1.0000    | —        | 0/2 ✗     | 0.5117      | 0.9065     | 0.8498       |
| 2%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0010     | 0.2406    | WSINDy   | 0/0 ✓     | 0.0016      | 0.2405     | 0.2399       |
| 2%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0007     | 0.2405    | WSINDy   | 0/0 ✓     | 0.0016      | 0.2405     | 0.2404       |
| 2%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0010     | 0.2406    | WSINDy   | 0/0 ✓     | 0.0016      | 0.2405     | 0.2399       |
| 2%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0007     | 0.2405    | WSINDy   | 0/0 ✓     | 0.0016      | 0.2405     | 0.2404       |
| 5%    | pointwise     | spectral   | —        | 0/4 ✗             | 1.0000     | 1.0000    | —        | 0/4 ✗     | 1.0000      | 1.0000     | 0.9699       |
| 5%    | pointwise     | polynomial | —        | 0/4 ✗             | 1.0000     | 1.0000    | —        | 0/4 ✗     | 1.0000      | 1.0000     | 0.9706       |
| 5%    | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0022     | 0.5229    | WSINDy   | 0/0 ✓     | 0.0043      | 0.5225     | 0.5211       |
| 5%    | weak (no IBP) | polynomial | 🟩        | 0/0 ✓             | 0.0013     | 0.5225    | WSINDy   | 0/0 ✓     | 0.0043      | 0.5225     | 0.5221       |
| 5%    | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0022     | 0.5229    | WSINDy   | 0/0 ✓     | 0.0043      | 0.5225     | 0.5211       |
| 5%    | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0013     | 0.5225    | WSINDy   | 0/0 ✓     | 0.0043      | 0.5225     | 0.5221       |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm   | d_tol/α   | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 2/4   | 2        | 0.2585  | 0.2565 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge     | 1.0       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.2       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.0005    | 0/0   | 2        | 0.0011  | 0.0019 | 0.1      | 🟢      |
| WSINDy      | Clean | RFE         |           | 0/0   | 2        | 0.0010  | 0.0019 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 3.0       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 2/4   | 2        | 0.9314  | 0.8647 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 1.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.005     | 0/0   | 4        | 0.0011  | 0.0019 | 0.3      | 🟢      |
| PolyJet     | Clean | RFE         |           | 0/0   | 4        | 0.0011  | 0.0019 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.3      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.005     | 0/0   | 2        | 0.0011  | 0.0019 | 0.0      | 🟢      |
| Multigraphs | Clean | RFE         |           | 0/0   | 2        | 0.0011  | 0.0019 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 0/0   | 2        | 0.0011  | 0.0019 | 19.3     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.0      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| Multigraphs | 2%    | Brute-force | L₀=1.0    | 1/4   | 1        | 0.9623  | 0.9306 | 19.0     | 🔴      |
| Multigraphs | 5%    | STRidge     | 3.0       | 0/4   | 0        | 1.0000  | 1.0000 | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| Multigraphs | 5%    | Brute-force | L₀=1.0    | 1/4   | 1        | 0.9901  | 0.9854 | 18.1     | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm   | d_tol/α   | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 2/4   | 2        | 0.2585  | 0.2565 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.0005    | 0/0   | 2        | 0.0013  | 0.0057 | 0.1      | 🟢      |
| WSINDy      | Clean | RFE         |           | 0/0   | 2        | 0.0012  | 0.0057 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.2       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| WSINDy      | 2%    | RFE         |           | 0/2   | 1        | 0.5117  | 0.9065 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 3.0       | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.01      | 0/0   | 4        | 0.0009  | 0.0051 | 0.2      | 🟢      |
| PolyJet     | Clean | RFE         |           | 0/0   | 4        | 0.0008  | 0.0051 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE         |           | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 1.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.005     | 0/0   | 2        | 0.0009  | 0.0051 | 0.0      | 🟢      |
| Multigraphs | Clean | RFE         |           | 0/0   | 2        | 0.0008  | 0.0051 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 0/0   | 2        | 0.0008  | 0.0051 | 19.2     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 0/4   | 0        | 1.0000  | 1.0000 | 0.0      | 🔴      |
| Multigraphs | 2%    | RFE         |           | 2/4   | 2        | 0.9285  | 0.8606 | 0.0      | 🔴      |
| Multigraphs | 2%    | Brute-force | L₀=1.0    | 0/2   | 1        | 0.4561  | 0.9244 | 18.5     | 🔴      |
| Multigraphs | 5%    | STRidge     | 0.2       | 0/4   | 0        | 1.0000  | 1.0000 | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE         |           | 1/4   | 1        | 0.9901  | 0.9854 | 0.0      | 🔴      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 1/2   | 2        | 0.9810  | 0.9722 | 18.1     | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm   | d_tol/α   | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 2.0       | 4  | no  | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | no  | 2/2   | 8        | 1.1109  | 0.1749 | 0.4      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 4  | no  | 2/4   | 2        | 0.2585  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | Clean | STRidge     | 0.05      | 4  | yes | 0/2   | 4        | 2.6929  | 0.5207 | 0.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | yes | 2/2   | 8        | 1.1109  | 0.1749 | 0.6      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 4  | yes | 2/4   | 2        | 0.2585  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge     | 0.05      | 4  | no  | 0/2   | 4        | 2.6112  | 0.7040 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | no  | 1/4   | 2        | 29.4599 | 0.4893 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 4  | no  | 2/4   | 2        | 0.9286  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge     | 0.05      | 4  | yes | 0/2   | 4        | 2.6112  | 0.7040 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | yes | 1/4   | 2        | 29.4599 | 0.4893 | 0.5      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 4  | yes | 2/4   | 2        | 0.9286  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 2.0       | 4  | no  | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 8  | no  | 2/2   | 6        | 0.6894  | 0.6660 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 4  | no  | 1/4   | 1        | 0.9902  | 0.9255 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 2.0       | 4  | yes | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 8  | yes | 2/2   | 6        | 0.6894  | 0.6660 | 0.5      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 4  | yes | 1/4   | 1        | 0.9902  | 0.9255 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.05      | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.05      | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.8360  | 0.8347 | 0.1      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.1       | 12 | no  | 0/0   | 2        | 0.0017  | 0.2406 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 4  | no  | 0/0   | 2        | 0.5144  | 0.7166 | 0.4      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | no  | 0/0   | 2        | 0.0016  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.1       | 12 | yes | 0/0   | 2        | 0.0017  | 0.2406 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 4  | yes | 0/0   | 2        | 0.5144  | 0.7166 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | yes | 0/0   | 2        | 0.0016  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 0.5       | 12 | no  | 0/0   | 2        | 0.0043  | 0.5225 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 8  | no  | 0/0   | 2        | 0.4008  | 0.7141 | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 12 | no  | 0/2   | 1        | 0.4924  | 0.7166 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.2       | 12 | yes | 0/0   | 2        | 0.0043  | 0.5225 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 8  | yes | 0/0   | 2        | 0.4008  | 0.7141 | 0.3      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 12 | yes | 0/2   | 1        | 0.4924  | 0.7166 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.5       | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.5      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge     | 0.5       | 12 | yes | 0/0   | 4        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0008  | 0.0006 | 0.4      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.2       | 8  | no  | 0/0   | 4        | 0.0310  | 0.3304 | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 8  | no  | 0/0   | 4        | 0.0023  | 0.3303 | 1.1      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 8  | no  | 0/1   | 3        | 1.2472  | 0.4159 | 0.0      | 🔴      |
| PolyJet     | 2%    | STRidge     | 0.2       | 8  | yes | 0/0   | 4        | 0.0310  | 0.3304 | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 8  | yes | 0/0   | 4        | 0.0023  | 0.3303 | 1.2      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 8  | yes | 0/1   | 3        | 1.2472  | 0.4159 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.5       | 8  | no  | 0/0   | 4        | 0.0673  | 0.6477 | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 8  | no  | 0/0   | 4        | 0.0045  | 0.6475 | 1.1      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 4  | no  | 1/4   | 1        | 0.9902  | 0.9255 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.5       | 8  | yes | 0/0   | 4        | 0.0673  | 0.6477 | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 8  | yes | 0/0   | 4        | 0.0045  | 0.6475 | 1.0      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 4  | yes | 1/4   | 1        | 0.9902  | 0.9255 | 0.0      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.5       | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 2        | 0.0008  | 0.0006 | 19.0     | 🟢      |
| Multigraphs | Clean | STRidge     | 0.5       | 12 | yes | 0/0   | 2        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.0008  | 0.0006 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | yes | 0/0   | 2        | 0.0008  | 0.0006 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 2        | 0.0008  | 0.0006 | 20.7     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 2        | 0.0010  | 0.2406 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.0020  | 0.2406 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | no  | 0/0   | 2        | 0.0020  | 0.2406 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 2        | 0.0020  | 0.2406 | 21.8     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 12 | yes | 0/0   | 2        | 0.0010  | 0.2406 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.0020  | 0.2406 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | yes | 0/0   | 2        | 0.0020  | 0.2406 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 2        | 0.0020  | 0.2406 | 18.2     | 🟢      |
| Multigraphs | 5%    | STRidge     | 1.0       | 8  | no  | 0/0   | 2        | 0.0032  | 0.6476 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.0022  | 0.5229 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 12 | no  | 0/0   | 2        | 0.0022  | 0.5229 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 12 | no  | 0/0   | 2        | 0.0022  | 0.5229 | 20.5     | 🟢      |
| Multigraphs | 5%    | STRidge     | 1.0       | 8  | yes | 0/0   | 2        | 0.0032  | 0.6476 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.0022  | 0.5229 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 12 | yes | 0/0   | 2        | 0.0022  | 0.5229 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 12 | yes | 0/0   | 2        | 0.0022  | 0.5229 | 18.9     | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm   | d_tol/α   | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 2.0       | 4  | no  | 0/4   | 0        | 1.0000  | 1.0000 | 0.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | no  | 2/2   | 8        | 1.1228  | 0.1769 | 0.4      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 4  | no  | 2/4   | 2        | 0.2585  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | Clean | STRidge     | 0.05      | 4  | yes | 0/2   | 4        | 2.6918  | 0.5208 | 0.1      | 🔴      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 4  | yes | 2/2   | 8        | 1.1228  | 0.1769 | 0.6      | 🔴      |
| PDE-FIND    | Clean | RFE         |           | 4  | yes | 2/4   | 2        | 0.2585  | 0.2518 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge     | 0.05      | 4  | no  | 0/2   | 4        | 2.5324  | 0.7057 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | no  | 1/4   | 3        | 44.1357 | 0.5277 | 0.4      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 4  | no  | 2/4   | 2        | 0.9286  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND    | 2%    | STRidge     | 0.05      | 4  | yes | 0/2   | 4        | 2.5324  | 0.7057 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | yes | 1/4   | 3        | 44.1357 | 0.5277 | 0.5      | 🔴      |
| PDE-FIND    | 2%    | RFE         |           | 4  | yes | 2/4   | 2        | 0.9286  | 0.5949 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.05      | 4  | no  | 0/2   | 4        | 2.4777  | 0.8937 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 4  | no  | 2/2   | 6        | 0.2303  | 0.8569 | 1.0      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 4  | no  | 1/4   | 1        | 0.9902  | 0.9255 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.05      | 4  | yes | 0/2   | 4        | 2.4777  | 0.8937 | 0.2      | 🔴      |
| PDE-FIND    | 5%    | MSTLS       | 0.75      | 4  | yes | 2/2   | 6        | 0.2303  | 0.8569 | 0.7      | 🔴      |
| PDE-FIND    | 5%    | RFE         |           | 4  | yes | 1/4   | 1        | 0.9902  | 0.9255 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.005     | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.8360  | 0.8347 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.005     | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.8360  | 0.8347 | 0.1      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.1       | 12 | no  | 0/0   | 2        | 0.0017  | 0.2406 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 4  | no  | 0/0   | 2        | 0.5209  | 0.7119 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | no  | 0/0   | 2        | 0.0016  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.1       | 12 | yes | 0/0   | 2        | 0.0017  | 0.2406 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 4  | yes | 0/0   | 2        | 0.5209  | 0.7119 | 0.3      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | yes | 0/0   | 2        | 0.0016  | 0.2405 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 0.2       | 12 | no  | 0/0   | 2        | 0.0043  | 0.5225 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 4  | no  | 0/0   | 2        | 0.0441  | 0.8533 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 12 | no  | 0/2   | 1        | 0.4924  | 0.7166 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.5       | 12 | yes | 0/0   | 2        | 0.0043  | 0.5225 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 4  | yes | 0/0   | 2        | 0.0441  | 0.8533 | 0.4      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 12 | yes | 0/2   | 1        | 0.4924  | 0.7166 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.5       | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.5      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge     | 0.5       | 12 | yes | 0/0   | 4        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0007  | 0.0007 | 0.5      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.5       | 8  | no  | 0/0   | 4        | 0.0304  | 0.3303 | 0.1      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 8  | no  | 0/0   | 4        | 0.0027  | 0.3302 | 0.9      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 8  | no  | 0/2   | 2        | 0.4695  | 0.4647 | 0.0      | 🔴      |
| PolyJet     | 2%    | STRidge     | 0.2       | 8  | yes | 0/0   | 4        | 0.0304  | 0.3303 | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 8  | yes | 0/0   | 4        | 0.0027  | 0.3302 | 0.9      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 8  | yes | 0/1   | 3        | 1.2410  | 0.4153 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.2       | 8  | no  | 0/0   | 4        | 0.0714  | 0.6471 | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 8  | no  | 0/0   | 4        | 0.0080  | 0.6469 | 0.9      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 8  | no  | 0/3   | 1        | 0.5488  | 0.8245 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.5       | 8  | yes | 0/0   | 4        | 0.0714  | 0.6471 | 0.2      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 8  | yes | 0/0   | 4        | 0.0080  | 0.6469 | 1.3      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 8  | yes | 0/3   | 1        | 0.5488  | 0.8245 | 0.0      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.1       | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 2        | 0.0007  | 0.0007 | 18.1     | 🟢      |
| Multigraphs | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 2        | 0.0007  | 0.0007 | 21.7     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 2        | 0.0007  | 0.2405 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 2        | 0.0024  | 0.2405 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | no  | 0/0   | 2        | 0.0024  | 0.2405 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 2        | 0.0024  | 0.2405 | 21.8     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 12 | yes | 0/0   | 2        | 0.0007  | 0.2405 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 2        | 0.0024  | 0.2405 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | yes | 0/0   | 2        | 0.0024  | 0.2405 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 2        | 0.0024  | 0.2405 | 20.1     | 🟢      |
| Multigraphs | 5%    | STRidge     | 2.0       | 12 | no  | 0/0   | 2        | 0.0013  | 0.5225 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 8  | no  | 0/0   | 2        | 0.0055  | 0.6469 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 8  | no  | 0/0   | 2        | 0.0055  | 0.6469 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 8  | no  | 0/0   | 2        | 0.0055  | 0.6469 | 22.1     | 🟢      |
| Multigraphs | 5%    | STRidge     | 2.0       | 12 | yes | 0/0   | 2        | 0.0013  | 0.5225 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 8  | yes | 0/0   | 2        | 0.0055  | 0.6469 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 8  | yes | 0/0   | 2        | 0.0055  | 0.6469 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.01   | 8  | yes | 0/0   | 2        | 0.0055  | 0.6469 | 20.7     | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)   |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:----------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0008       | 1144526.1 |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0003       | 1933914.7 |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0000       | 3063426.0 |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0012       | 399891.4  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0003       | 1103259.0 |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0001       | 1896428.8 |
| Clean | weak (IBP)    | spectral   | 4  | yes | 0.0008       | 3534221.1 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0003       | 1940044.0 |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0000       | 3065203.8 |
| Clean | weak (IBP)    | polynomial | 4  | yes | 0.0013       | 1438261.1 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0003       | 1060615.2 |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0001       | 1884770.8 |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.5581       | 5461.2    |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.3300       | 14793.5   |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.2399       | 25197.9   |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.5568       | 26996.4   |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.3300       | 74434.3   |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.2404       | 174841.5  |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 0.5581       | 18061.6   |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.3300       | 14794.1   |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.2399       | 25197.6   |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 0.5568       | 88760.8   |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.3300       | 72730.0   |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.2404       | 174824.0  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.8498       | 2218.6    |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.6465       | 5152.9    |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.5211       | 11056.3   |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.8488       | 9053.8    |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.6466       | 30716.2   |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.5221       | 62737.0   |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 0.8499       | 4890.3    |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.6465       | 5133.2    |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.5211       | 11050.3   |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 0.8487       | 31368.5   |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.6466       | 30454.4   |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.5221       | 62734.9   |
---
| Source                                      | Wall time | CPU                                      | Threads |
|:--------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260710_003218_pme_mstls_fixed.csv` | 97.1 min  | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 4       |
