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


| Noise | Method        | Backend    | Eq wins? | Multigraphs fp/fn | Multigraphs E_coeff | Mult. E_data | Non-eq best | Non-eq fp/fn | Non-eq E_coeff | Non-eq E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:--------------------|:-------------|:------------|:-------------|:---------------|:--------------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0080              | 0.0664       | PolyJet     | 0/0 ✓        | 0.0080         | 0.0664        | 0.0441       |
| Clean | pointwise     | polynomial | ⚪        | 0/0 ✓             | 0.0114              | 0.1133       | PolyJet     | 0/0 ✓        | 0.0115         | 0.1133        | 0.1075       |
| Clean | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0048              | 0.0032       | PolyJet     | 0/0 ✓        | 0.0048         | 0.0032        | 0.0028       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0048              | 0.0032       | PolyJet     | 0/0 ✓        | 0.0048         | 0.0032        | 0.0028       |
| Clean | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0048              | 0.0032       | PDE-FIND    | 0/0 ✓        | 0.0048         | 0.0032        | 0.0028       |
| Clean | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0048              | 0.0032       | PDE-FIND    | 0/0 ✓        | 0.0048         | 0.0032        | 0.0028       |
| 2%    | pointwise     | spectral   | —        | 0/2 ✗             | 0.7637              | 0.4072       | —           | 0/2 ✗        | 0.7637         | 0.4072        | 0.2885       |
| 2%    | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0515              | 0.2545       | PDE-FIND    | 0/0 ✓        | 0.0521         | 0.2545        | 0.2412       |
| 2%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0050              | 0.0217       | PolyJet     | 0/0 ✓        | 0.0050         | 0.0217        | 0.0216       |
| 2%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0050              | 0.0217       | PolyJet     | 0/0 ✓        | 0.0049         | 0.0217        | 0.0216       |
| 2%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0050              | 0.0217       | PolyJet     | 0/0 ✓        | 0.0050         | 0.0217        | 0.0216       |
| 2%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0050              | 0.0217       | PolyJet     | 0/0 ✓        | 0.0050         | 0.0217        | 0.0216       |
| 5%    | pointwise     | spectral   | —        | 0/2 ✗             | 0.9539              | 0.5669       | —           | 0/2 ✗        | 0.9539         | 0.5669        | 0.5132       |
| 5%    | pointwise     | polynomial | 🔴        | 0/2 ✗             | 0.9539              | 0.5669       | PDE-FIND    | 0/0 ✓        | 0.3388         | 0.5149        | 0.4921       |
| 5%    | weak (no IBP) | spectral   | ⚪        | 0/0 ✓             | 0.0070              | 0.0536       | PolyJet     | 0/0 ✓        | 0.0069         | 0.0536        | 0.0534       |
| 5%    | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0069              | 0.0536       | PDE-FIND    | 0/0 ✓        | 0.0068         | 0.0536        | 0.0534       |
| 5%    | weak (IBP)    | spectral   | ⚪        | 0/0 ✓             | 0.0070              | 0.0536       | PDE-FIND    | 0/0 ✓        | 0.0069         | 0.0536        | 0.0534       |
| 5%    | weak (IBP)    | polynomial | ⚪        | 0/0 ✓             | 0.0070              | 0.0536       | PDE-FIND    | 0/0 ✓        | 0.0069         | 0.0536        | 0.0534       |

## Pointwise formulation

### Spectral derivatives

| Library     | Noise | Algorithm   | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.05    | 0/0   | 4        | 0.0084  | 0.0664 | 0.1      | 🟢      |
| PDE-FIND    | Clean | RFE         |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.2     | 0/2   | 2        | 0.7637  | 0.4072 | 0.1      | 🔴      |
| PDE-FIND    | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.2     | 0/2   | 2        | 0.9539  | 0.5669 | 0.1      | 🔴      |
| PDE-FIND    | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.05    | 0/0   | 4        | 0.0084  | 0.0664 | 0.1      | 🟢      |
| WSINDy      | Clean | RFE         |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.2     | 0/2   | 2        | 0.7637  | 0.4072 | 0.1      | 🔴      |
| WSINDy      | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.2     | 0/2   | 2        | 0.9539  | 0.5669 | 0.1      | 🔴      |
| WSINDy      | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.05    | 0/0   | 4        | 0.0084  | 0.0664 | 0.2      | 🟢      |
| PolyJet     | Clean | RFE         |         | 0/0   | 4        | 0.0080  | 0.0664 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.2     | 0/2   | 2        | 0.7637  | 0.4072 | 0.2      | 🔴      |
| PolyJet     | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.2     | 0/2   | 2        | 0.9539  | 0.5669 | 0.2      | 🔴      |
| PolyJet     | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.05    | 0/0   | 3        | 0.0084  | 0.0664 | 0.0      | 🟢      |
| Multigraphs | Clean | RFE         |         | 0/0   | 3        | 0.0080  | 0.0664 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.01 | 0/0   | 3        | 0.0080  | 0.0664 | 19.1     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.2     | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| Multigraphs | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| Multigraphs | 2%    | Brute-force | L₀=1.0  | 0/3   | 1        | 0.8719  | 0.7158 | 17.9     | 🔴      |
| Multigraphs | 5%    | STRidge     | 0.2     | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| Multigraphs | 5%    | Brute-force | L₀=1.0  | 0/3   | 1        | 0.9758  | 0.7732 | 18.1     | 🔴      |

### Polynomial stencils

| Library     | Noise | Algorithm   | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.05    | 0/0   | 4        | 0.0117  | 0.1133 | 0.1      | 🟢      |
| PDE-FIND    | Clean | RFE         |         | 0/0   | 4        | 0.0115  | 0.1133 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.05    | 0/0   | 4        | 0.0521  | 0.2545 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.05    | 0/0   | 4        | 0.3388  | 0.5149 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.05    | 0/0   | 4        | 0.0117  | 0.1133 | 0.1      | 🟢      |
| WSINDy      | Clean | RFE         |         | 0/0   | 4        | 0.0115  | 0.1133 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.05    | 0/0   | 4        | 0.0521  | 0.2545 | 0.1      | 🟢      |
| WSINDy      | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.05    | 0/0   | 4        | 0.3388  | 0.5149 | 0.1      | 🟢      |
| WSINDy      | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.01    | 0/0   | 4        | 0.0117  | 0.1133 | 0.2      | 🟢      |
| PolyJet     | Clean | RFE         |         | 0/0   | 4        | 0.0115  | 0.1133 | 0.1      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.05    | 0/0   | 4        | 0.0521  | 0.2545 | 0.2      | 🟢      |
| PolyJet     | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.05    | 0/0   | 4        | 0.3388  | 0.5149 | 0.2      | 🟢      |
| PolyJet     | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.05    | 0/0   | 3        | 0.0116  | 0.1133 | 0.0      | 🟢      |
| Multigraphs | Clean | RFE         |         | 0/0   | 3        | 0.0114  | 0.1133 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.01 | 0/0   | 3        | 0.0114  | 0.1133 | 19.2     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.05    | 0/0   | 3        | 0.0515  | 0.2545 | 0.0      | 🟢      |
| Multigraphs | 2%    | RFE         |         | 0/2   | 2        | 0.7637  | 0.4072 | 0.0      | 🔴      |
| Multigraphs | 2%    | Brute-force | L₀=0.01 | 0/0   | 3        | 0.0540  | 0.2545 | 18.0     | 🟢      |
| Multigraphs | 5%    | STRidge     | 0.2     | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| Multigraphs | 5%    | RFE         |         | 0/2   | 2        | 0.9539  | 0.5669 | 0.0      | 🔴      |
| Multigraphs | 5%    | Brute-force | L₀=1.0  | 0/3   | 1        | 0.9758  | 0.7732 | 18.0     | 🔴      |

## Weak formulation

### Spectral derivatives

| Library     | Noise | Algorithm   | d_tol/α   | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.05      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND    | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND    | Clean | STRidge     | 0.05      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND    | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.2       | 12 | no  | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND    | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND    | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge     | 0.2       | 12 | no  | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0078  | 0.0536 | 0.3      | 🟢      |
| PDE-FIND    | 5%    | RFE         |           | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0078  | 0.0536 | 0.4      | 🟢      |
| PDE-FIND    | 5%    | RFE         |           | 4  | yes | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.05      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 2.8796  | 2.8813 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.05      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 2.8796  | 2.8813 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.05      | 12 | no  | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 4.9910  | 2.8793 | 0.2      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.05      | 12 | yes | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 4.9910  | 2.8793 | 0.2      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 0.2       | 12 | no  | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 6.0213  | 2.8685 | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 6.0212  | 2.8685 | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 4  | yes | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.2       | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.1       | 12 | no  | 0/0   | 4        | 0.0050  | 0.0217 | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.1       | 12 | yes | 0/0   | 4        | 0.0050  | 0.0217 | 0.2      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PolyJet     | 5%    | STRidge     | 0.5       | 12 | no  | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0078  | 0.0536 | 0.9      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.5       | 12 | yes | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0078  | 0.0536 | 0.9      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 4  | yes | 0/2   | 2        | 0.9545  | 0.3345 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.1       | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 19.3     | 🟢      |
| Multigraphs | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 19.2     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.1       | 12 | no  | 0/0   | 3        | 0.0050  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 18.0     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.1       | 12 | yes | 0/0   | 3        | 0.0050  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 3        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | yes | 0/0   | 3        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 3        | 0.0052  | 0.0217 | 18.9     | 🟢      |
| Multigraphs | 5%    | STRidge     | 1.5       | 12 | no  | 0/0   | 3        | 0.0070  | 0.0536 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 3        | 0.0078  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 8  | no  | 0/0   | 3        | 0.0103  | 0.0770 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 3        | 0.0078  | 0.0536 | 19.0     | 🟢      |
| Multigraphs | 5%    | STRidge     | 1.5       | 12 | yes | 0/0   | 3        | 0.0070  | 0.0536 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 3        | 0.0078  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 8  | yes | 0/0   | 3        | 0.0104  | 0.0770 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 3        | 0.0078  | 0.0536 | 19.0     | 🟢      |

### Polynomial stencils

| Library     | Noise | Algorithm   | d_tol/α   | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:------------|:------|:------------|:----------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PDE-FIND    | Clean | STRidge     | 0.1       | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND    | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND    | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PDE-FIND    | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.3      | 🟢      |
| PDE-FIND    | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.1       | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0051  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND    | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0051  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND    | 2%    | STRidge     | 0.1       | 12 | yes | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| PDE-FIND    | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.3      | 🟢      |
| PDE-FIND    | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PDE-FIND    | 5%    | STRidge     | 0.2       | 12 | no  | 0/0   | 4        | 0.0068  | 0.0536 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0077  | 0.0536 | 0.3      | 🟢      |
| PDE-FIND    | 5%    | RFE         |           | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PDE-FIND    | 5%    | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PDE-FIND    | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0078  | 0.0536 | 0.3      | 🟢      |
| PDE-FIND    | 5%    | RFE         |           | 4  | yes | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| WSINDy      | Clean | STRidge     | 0.05      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 2.8797  | 2.8812 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy      | Clean | STRidge     | 0.05      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| WSINDy      | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 2.8796  | 2.8813 | 0.2      | 🟢      |
| WSINDy      | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.05      | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 4.9913  | 2.8792 | 0.2      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0051  | 0.0217 | 0.0      | 🟢      |
| WSINDy      | 2%    | STRidge     | 0.05      | 12 | yes | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| WSINDy      | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 4.9910  | 2.8793 | 0.2      | 🟢      |
| WSINDy      | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| WSINDy      | 5%    | STRidge     | 0.2       | 12 | no  | 0/0   | 4        | 0.0068  | 0.0536 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 6.0216  | 2.8685 | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| WSINDy      | 5%    | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| WSINDy      | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 6.0212  | 2.8685 | 0.2      | 🟢      |
| WSINDy      | 5%    | RFE         |           | 4  | yes | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PolyJet     | Clean | STRidge     | 0.05      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.2      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.7      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | no  | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet     | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| PolyJet     | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.8      | 🟢      |
| PolyJet     | Clean | RFE         |           | 12 | yes | 0/0   | 4        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 4        | 0.0049  | 0.0217 | 0.1      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0051  | 0.0217 | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 12 | no  | 0/0   | 4        | 0.0051  | 0.0217 | 0.0      | 🟢      |
| PolyJet     | 2%    | STRidge     | 0.2       | 12 | yes | 0/0   | 4        | 0.0050  | 0.0217 | 0.1      | 🟢      |
| PolyJet     | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.8      | 🟢      |
| PolyJet     | 2%    | RFE         |           | 12 | yes | 0/0   | 4        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| PolyJet     | 5%    | STRidge     | 0.5       | 12 | no  | 0/0   | 4        | 0.0068  | 0.0536 | 0.1      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 4        | 0.0077  | 0.0536 | 1.0      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 4  | no  | 0/2   | 2        | 0.9545  | 0.3345 | 0.0      | 🔴      |
| PolyJet     | 5%    | STRidge     | 0.5       | 12 | yes | 0/0   | 4        | 0.0069  | 0.0536 | 0.1      | 🟢      |
| PolyJet     | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 4        | 0.0078  | 0.0536 | 1.0      | 🟢      |
| PolyJet     | 5%    | RFE         |           | 4  | yes | 0/2   | 2        | 0.9545  | 0.3345 | 0.1      | 🔴      |
| Multigraphs | Clean | STRidge     | 0.2       | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 3        | 0.0048  | 0.0032 | 18.3     | 🟢      |
| Multigraphs | Clean | STRidge     | 0.1       | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | MSTLS       | 0.25      | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.1      | 🟢      |
| Multigraphs | Clean | RFE         |           | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 0.0      | 🟢      |
| Multigraphs | Clean | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 3        | 0.0048  | 0.0032 | 18.3     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.5       | 12 | no  | 0/0   | 3        | 0.0050  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 3        | 0.0052  | 0.0217 | 18.4     | 🟢      |
| Multigraphs | 2%    | STRidge     | 0.1       | 12 | yes | 0/0   | 3        | 0.0050  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 3        | 0.0052  | 0.0217 | 0.1      | 🟢      |
| Multigraphs | 2%    | RFE         |           | 12 | yes | 0/0   | 3        | 0.0052  | 0.0217 | 0.0      | 🟢      |
| Multigraphs | 2%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 3        | 0.0052  | 0.0217 | 19.2     | 🟢      |
| Multigraphs | 5%    | STRidge     | 0.5       | 12 | no  | 0/0   | 3        | 0.0069  | 0.0536 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | no  | 0/0   | 3        | 0.0077  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 8  | no  | 0/0   | 3        | 0.0096  | 0.0770 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.0001 | 12 | no  | 0/0   | 3        | 0.0077  | 0.0536 | 18.9     | 🟢      |
| Multigraphs | 5%    | STRidge     | 0.5       | 12 | yes | 0/0   | 3        | 0.0070  | 0.0536 | 0.0      | 🟢      |
| Multigraphs | 5%    | MSTLS       | 0.25      | 12 | yes | 0/0   | 3        | 0.0078  | 0.0536 | 0.1      | 🟢      |
| Multigraphs | 5%    | RFE         |           | 8  | yes | 0/0   | 3        | 0.0104  | 0.0770 | 0.0      | 🟢      |
| Multigraphs | 5%    | Brute-force | L₀=0.0001 | 12 | yes | 0/0   | 3        | 0.0078  | 0.0536 | 19.1     | 🟢      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Method        | Backend    | m  | IBP | E_data floor | cond(Θ) |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:--------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0035       | 1742.0  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0030       | 1475.2  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0028       | 2683.0  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0112       | 1704.3  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0030       | 1418.8  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0028       | 2568.6  |
| Clean | weak (IBP)    | spectral   | 4  | yes | 0.0035       | 25684.2 |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.0030       | 1477.8  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0028       | 2682.2  |
| Clean | weak (IBP)    | polynomial | 4  | yes | 0.0035       | 25885.5 |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.0029       | 1406.6  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0028       | 2563.8  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 0.0621       | 2330.4  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 0.0309       | 2821.4  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 0.0216       | 4037.5  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 0.0618       | 1731.3  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.0309       | 1427.5  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.0216       | 2582.2  |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 0.0625       | 11816.9 |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 0.0309       | 2841.3  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 0.0216       | 4037.9  |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 0.0626       | 11217.1 |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 0.0309       | 1414.7  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.0216       | 2577.5  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 0.1502       | 10382.8 |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 0.0766       | 13252.1 |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 0.0534       | 19563.0 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 0.1467       | 2046.5  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 0.0766       | 1462.6  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 0.0534       | 2531.0  |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 0.1515       | 12502.2 |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 0.0766       | 13321.0 |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 0.0534       | 19566.7 |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 0.1500       | 9239.9  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 0.0766       | 1449.7  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 0.0534       | 2527.3  |
---
| Source                                     | Wall time | CPU                                      | Threads |
|:-------------------------------------------|----------:|:-----------------------------------------|:-------:|
| `sweep_20260709_144235_ac_mstls_fixed.csv` | 100.6 min | Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz | 1       |
