# Kuramoto-Sivashinsky equation

The equation under consideration is

$$f_t = -\Delta f - \Delta^2 f - \tfrac12\|\nabla f\|^2$$

which can be written in terms of $O(2)$-equivariant terms

| Operator         | True coefficient |
|:-----------------|:-----------------|
| $\Delta f$       | $-1.0$           |
| $\Delta^2 f$     | $-1.0$           |
| $\|\nabla f\|^2$ | $-0.5$           |

or standard terms (PolyJet basis)

| Operator            | True coeff | PDE-FIND | WSINDy | Poly jet | Multigraph |
|:--------------------|:-----------|:---------|:-------|:---------|:-----------|
| $\partial_{xx} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yy} f$   |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxxx} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{xxyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $\partial_{yyyy} f$ |            | 🟢        | 🟢      | 🟢        | 🟢          |
| $f_x^2$             |            | 🔴        | 🔴      | 🟢        | 🟢          |
| $f_y^2$             |            | 🔴        | 🔴      | 🟢        | 🟢          |

KS contains the gradient norm $\|\nabla f\|^2 = f_x^2 + f_y^2$ which neither PDE-FIND nor WSINDy can represent as a sparse column (PDE-FIND lacks products of derivatives; in WSINDy, $f_x^2$ is locked with $f\cdot f_{xx}$ inside $\partial_{xx}(f^2)$).  Three KS regimes are simulated: cellular ($L=4\pi$), weak-chaotic ($L=8\pi$), chaotic ($L=16\pi$).

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


| Noise | Formulation   | Backend    | Eq wins? | Multigraphs fp/fn | Mg E_coeff | Mult. best             | Mg E_data | Std best      | Std fp/fn | Std E_coeff | Std E_data | E_data floor |
|:------|:-------------:|:----------:|:--------:|:------------------|:-----------|:-----------------------|:----------|:--------------|:----------|:------------|:-----------|-------------:|
| Clean | pointwise     | spectral   | ⚪        | 0/0 ✓             | 0.0047     | multigraph/RFE         | 0.0176    | PolyJet/RFE   | 0/0 ✓     | 0.0047      | 0.0175     | 0.0182       |
| Clean | pointwise     | polynomial | 🟩        | 0/0 ✓             | 0.0117     | multigraph/RFE         | 0.0325    | PolyJet/RFE   | 0/0 ✓     | 0.0128      | 0.0320     | 0.0345       |
| Clean | weak (no IBP) | spectral   | 🟩        | 0/0 ✓             | 0.0021     | mod_multigraph/STRidge | 0.0112    | PolyJet/MSTLS | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (no IBP) | polynomial | ⚪        | 0/0 ✓             | 0.0022     | multigraph/RFE         | 0.0112    | PolyJet/RFE   | 0/0 ✓     | 0.0023      | 0.0112     | 0.0114       |
| Clean | weak (IBP)    | spectral   | 🟩        | 0/0 ✓             | 0.0279     | multigraph/MSTLS       | 0.0158    | PolyJet/MSTLS | 0/0 ✓     | 0.0280      | 0.0157     | 0.0326       |
| Clean | weak (IBP)    | polynomial | 🟩        | 0/0 ✓             | 0.0278     | multigraph/RFE         | 0.0158    | PolyJet/MSTLS | 0/0 ✓     | 0.0279      | 0.0157     | 0.0326       |
| 2%    | pointwise     | spectral   | —        | 18/0 ✗            | 1.0000     | —                      | 0.6378    | —             | 73/0 ✗    | 1.0000      | 0.6364     | 1381.7721    |
| 2%    | pointwise     | polynomial | —        | 29/0 ✗            | 0.9994     | —                      | 0.6238    | —             | 97/0 ✗    | 0.9991      | 0.6182     | 29.8498      |
| 2%    | weak (no IBP) | spectral   | —        | 6/0 ✗             | 0.9307     | —                      | 0.3476    | —             | 66/0 ✗    | 0.7213      | 0.2564     | 11.3180      |
| 2%    | weak (no IBP) | polynomial | —        | 1/0 ✗             | 0.3372     | —                      | 0.1866    | —             | 9/0 ✗     | 0.1883      | 0.1099     | 0.8842       |
| 2%    | weak (IBP)    | spectral   | —        | 9/0 ✗             | 0.8618     | —                      | 0.3082    | —             | 11/0 ✗    | 0.9200      | 0.3359     | 11.4293      |
| 2%    | weak (IBP)    | polynomial | —        | 1/0 ✗             | 0.3329     | —                      | 0.1858    | —             | 9/0 ✗     | 0.1972      | 0.1098     | 0.9103       |
| 5%    | pointwise     | spectral   | —        | 34/0 ✗            | 1.0000     | —                      | 0.8725    | —             | 86/0 ✗    | 1.0000      | 0.8714     | 2190.0271    |
| 5%    | pointwise     | polynomial | —        | 31/0 ✗            | 0.9996     | —                      | 0.8717    | —             | 86/0 ✗    | 0.9994      | 0.8702     | 47.4252      |
| 5%    | weak (no IBP) | spectral   | —        | 28/0 ✗            | 1.0001     | —                      | 0.5369    | —             | 42/0 ✗    | 0.9995      | 0.4133     | 69.4094      |
| 5%    | weak (no IBP) | polynomial | —        | 5/0 ✗             | 0.9730     | —                      | 0.3915    | —             | 13/0 ✗    | 0.6399      | 0.2611     | 5.4493       |
| 5%    | weak (IBP)    | spectral   | —        | 24/0 ✗            | 1.0001     | —                      | 0.5370    | —             | 58/0 ✗    | 1.0000      | 0.5348     | 70.0482      |
| 5%    | weak (IBP)    | polynomial | —        | 5/0 ✗             | 0.9617     | —                      | 0.3901    | —             | 13/0 ✗    | 0.6181      | 0.2604     | 5.4753       |


## Ablation

Each row shows the best configuration for a (library, algorithm) pair, ranked from best to worst (exact recovery, then fn, then fp, then E_coeff).  vs. best = (E − E_best) / E_best × 100% (✗ if method does not achieve exact recovery).

### Clean (0% noise)

| Method                         | Best conf.                    | fp/fn | E_coeff | Time (s) | vs. best |
|:-------------------------------|:------------------------------|:------|--------:|---------:|---------:|
| Modified Multigraphs + STRidge | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0021  | 0.2      | —        |
| Multigraphs + MSTLS            | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0022  | 0.3      | +3.6%    |
| Modified Multigraphs + MSTLS   | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0022  | 0.3      | +3.6%    |
| Multigraphs + RFE              | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0022  | 0.2      | +3.6%    |
| Modified Multigraphs + RFE     | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0022  | 0.2      | +3.6%    |
| PolyJet + MSTLS                | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0023  | 2.4      | +8.3%    |
| PolyJet + RFE                  | spectral, weak (no IBP), m=12 | 0/0 ✓ | 0.0023  | 0.2      | +8.3%    |
| Multigraphs + STRidge          | spectral, weak (no IBP), m=8  | 0/0 ✓ | 0.0031  | 0.2      | +43.1%   |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8  | 0/0 ✓ | 0.0032  | 0.8      | +48.6%   |

### 2% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + MSTLS            | polynomial, weak (IBP), m=12    | 1/0 ✗  | 0.3329  | 0.3      | ✗        |
| Modified Multigraphs + MSTLS   | polynomial, weak (IBP), m=12    | 1/0 ✗  | 0.3329  | 0.3      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (no IBP), m=12 | 9/0 ✗  | 0.1883  | 2.2      | ✗        |
| PolyJet + STRidge              | polynomial, weak (IBP), m=12    | 16/0 ✗ | 0.2185  | 0.7      | ✗        |
| Multigraphs + STRidge          | spectral, pointwise             | 18/0 ✗ | 1.0000  | 0.3      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (no IBP), m=12 | 27/0 ✗ | 0.2379  | 0.2      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9984  | 0.2      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9984  | 0.2      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=8  | 0/6 ✗  | 0.9987  | 0.4      | ✗        |

### 5% noise

| Method                         | Best conf.                      | fp/fn  | E_coeff | Time (s) | vs. best |
|:-------------------------------|:--------------------------------|:-------|--------:|---------:|---------:|
| Multigraphs + MSTLS            | polynomial, weak (IBP), m=8     | 5/0 ✗  | 0.9617  | 0.4      | ✗        |
| Modified Multigraphs + MSTLS   | polynomial, weak (no IBP), m=12 | 10/0 ✗ | 0.8007  | 0.3      | ✗        |
| PolyJet + MSTLS                | polynomial, weak (IBP), m=12    | 13/0 ✗ | 0.6181  | 2.3      | ✗        |
| Multigraphs + STRidge          | spectral, weak (IBP), m=4       | 24/0 ✗ | 1.0001  | 0.3      | ✗        |
| Modified Multigraphs + STRidge | polynomial, weak (IBP), m=12    | 30/0 ✗ | 0.9072  | 0.2      | ✗        |
| PolyJet + STRidge              | spectral, weak (no IBP), m=8    | 42/0 ✗ | 0.9995  | 1.1      | ✗        |
| Modified Multigraphs + RFE     | polynomial, weak (no IBP), m=12 | 0/5 ✗  | 0.9988  | 0.2      | ✗        |
| Multigraphs + RFE              | polynomial, weak (no IBP), m=8  | 0/5 ✗  | 0.9988  | 0.2      | ✗        |
| PolyJet + RFE                  | polynomial, weak (no IBP), m=12 | 1/7 ✗  | 0.9996  | 0.3      | ✗        |

## Pointwise formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0047  | 0.0175 | 1.2      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0047  | 0.0175 | 0.6      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.1     | 73/0  | 80       | 1.0000  | 0.6364 | 1.7      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.5      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 86/0  | 93       | 1.0000  | 0.8714 | 1.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.5      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 0/0   | 3        | 0.0047  | 0.0176 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 18/0  | 12       | 1.0000  | 0.6378 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.1     | 34/0  | 18       | 1.0000  | 0.8725 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 0/0   | 3        | 0.0047  | 0.0176 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0047  | 0.0176 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.01    | 35/0  | 18       | 1.0000  | 0.6377 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.6433 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 1.0000  | 0.8725 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 1.0000  | 0.8979 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 1.0     | 0/0   | 7        | 0.0130  | 0.0320 | 1.0      | 🟢      |
| PolyJet          | Clean | RFE       |         | 0/0   | 7        | 0.0128  | 0.0320 | 0.5      | 🟢      |
| PolyJet          | 2%    | STRidge   | 0.001   | 97/0  | 104      | 0.9991  | 0.6182 | 1.9      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.5      | 🔴      |
| PolyJet          | 5%    | STRidge   | 0.1     | 86/0  | 93       | 0.9994  | 0.8702 | 1.9      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.5      | 🔴      |
| Multigraphs      | Clean | STRidge   | 1.0     | 0/0   | 3        | 0.0119  | 0.0325 | 0.2      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 0/0   | 3        | 0.0117  | 0.0325 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 0.1     | 29/0  | 16       | 0.9994  | 0.6238 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.01    | 31/0  | 17       | 0.9996  | 0.8717 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 0/0   | 3        | 0.0119  | 0.0325 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 0/0   | 3        | 0.0117  | 0.0325 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 0.1     | 31/0  | 17       | 0.9994  | 0.6238 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 1/7   | 1        | 0.9991  | 0.6433 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.001   | 37/0  | 20       | 0.9996  | 0.8717 | 0.4      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 1/7   | 1        | 0.9996  | 0.8979 | 0.2      | 🔴      |

## Weak formulation

### Spectral derivatives

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 7        | 0.0032  | 0.0138 | 0.8      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023  | 0.0112 | 2.4      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023  | 0.0112 | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 8        | 0.0297  | 0.0100 | 0.6      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0280  | 0.0157 | 2.1      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0280  | 0.0157 | 0.3      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | no  | 66/0  | 73       | 0.7213  | 0.2564 | 0.8      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.75    | 8  | no  | 12/1  | 18       | 0.9369  | 0.3362 | 3.0      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | no  | 0/6   | 1        | 1.0000  | 0.3548 | 0.2      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 66/0  | 73       | 0.7214  | 0.2564 | 1.0      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 8  | yes | 11/0  | 18       | 0.9200  | 0.3359 | 2.7      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 12 | yes | 0/6   | 1        | 1.0000  | 0.3548 | 0.2      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 8  | no  | 42/0  | 49       | 0.9995  | 0.4133 | 1.1      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | no  | 9/3   | 13       | 1.0009  | 0.3575 | 2.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000  | 0.6460 | 0.5      | 🔴      |
| PolyJet          | 5%    | STRidge   | 1.0     | 4  | yes | 58/0  | 65       | 1.0000  | 0.5348 | 1.7      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.75    | 12 | yes | 9/3   | 13       | 1.0008  | 0.3575 | 2.1      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 4  | yes | 1/7   | 1        | 1.0000  | 0.6460 | 0.5      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0031  | 0.0139 | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.2      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.0297  | 0.0100 | 0.2      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0279  | 0.0158 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0279  | 0.0158 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 32/0  | 18       | 0.7870  | 0.2829 | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 8  | no  | 6/0   | 6        | 0.9307  | 0.3476 | 0.4      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 1.0000  | 0.3445 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 4  | yes | 27/0  | 16       | 1.0002  | 0.4177 | 0.3      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.75    | 12 | yes | 14/0  | 9        | 0.7357  | 0.2872 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 1.0000  | 0.3445 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | no  | 28/0  | 17       | 1.0001  | 0.5369 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 1.0000  | 0.3701 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000  | 0.6460 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | yes | 24/0  | 14       | 1.0001  | 0.5370 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.25    | 12 | yes | 1/7   | 1        | 1.0000  | 0.3701 | 0.3      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 4  | yes | 1/7   | 1        | 1.0000  | 0.6460 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0021  | 0.0112 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.0297  | 0.0100 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0279  | 0.0158 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0279  | 0.0158 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 12 | no  | 34/0  | 19       | 0.7852  | 0.2827 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | no  | 9/0   | 7        | 0.7235  | 0.2876 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 1.0000  | 0.3445 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 1.0     | 4  | yes | 31/0  | 15       | 1.0001  | 0.4268 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.75    | 12 | yes | 9/0   | 7        | 0.8618  | 0.3082 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 1.0000  | 0.3445 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | no  | 35/0  | 18       | 1.0295  | 0.3541 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | no  | 1/7   | 1        | 1.0000  | 0.3701 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | no  | 1/7   | 1        | 1.0000  | 0.6460 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 10.0    | 12 | yes | 31/0  | 14       | 1.0336  | 0.3544 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 4  | yes | 3/3   | 5        | 1.0000  | 0.5377 | 0.5      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 4  | yes | 1/7   | 1        | 1.0000  | 0.6460 | 0.2      | 🔴      |

### Polynomial stencils

| Library          | Noise | Algorithm | d_tol/α | m  | IBP | fp/fn | n_active | E_coeff | E_data | Time (s) | Exact? |
|:-----------------|:------|:----------|:--------|:---|:----|:------|:---------|:--------|:-------|:---------|:-------|
| PolyJet          | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 7        | 0.0042  | 0.0141 | 0.8      | 🟢      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | no  | 0/0   | 7        | 0.0023  | 0.0112 | 1.9      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | no  | 0/0   | 7        | 0.0023  | 0.0112 | 0.2      | 🟢      |
| PolyJet          | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 8        | 0.0297  | 0.0100 | 0.6      | 🔴      |
| PolyJet          | Clean | MSTLS     | 0.25    | 12 | yes | 0/0   | 7        | 0.0279  | 0.0157 | 2.3      | 🟢      |
| PolyJet          | Clean | RFE       |         | 12 | yes | 0/0   | 7        | 0.0279  | 0.0157 | 0.3      | 🟢      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | no  | 17/0  | 24       | 0.2182  | 0.1089 | 0.8      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | no  | 9/0   | 16       | 0.1883  | 0.1099 | 2.2      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | no  | 0/6   | 1        | 0.9987  | 0.4262 | 0.4      | 🔴      |
| PolyJet          | 2%    | STRidge   | 10.0    | 12 | yes | 16/0  | 23       | 0.2185  | 0.1091 | 0.7      | 🔴      |
| PolyJet          | 2%    | MSTLS     | 0.5     | 12 | yes | 9/0   | 16       | 0.1972  | 0.1098 | 2.1      | 🔴      |
| PolyJet          | 2%    | RFE       |         | 8  | yes | 0/6   | 1        | 0.9987  | 0.4262 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | no  | 55/0  | 62       | 0.6891  | 0.2526 | 0.9      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | no  | 13/0  | 20       | 0.6399  | 0.2611 | 2.2      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | no  | 1/7   | 1        | 0.9996  | 0.5246 | 0.3      | 🔴      |
| PolyJet          | 5%    | STRidge   | 10.0    | 12 | yes | 54/0  | 61       | 0.6911  | 0.2526 | 0.9      | 🔴      |
| PolyJet          | 5%    | MSTLS     | 0.5     | 12 | yes | 13/0  | 20       | 0.6181  | 0.2604 | 2.3      | 🔴      |
| PolyJet          | 5%    | RFE       |         | 12 | yes | 1/7   | 1        | 0.9996  | 0.5246 | 0.3      | 🔴      |
| Multigraphs      | Clean | STRidge   | 10.0    | 8  | no  | 0/0   | 3        | 0.0039  | 0.0141 | 0.2      | 🟢      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.2      | 🟢      |
| Multigraphs      | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.0296  | 0.0100 | 0.2      | 🔴      |
| Multigraphs      | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0278  | 0.0158 | 0.3      | 🟢      |
| Multigraphs      | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0278  | 0.0158 | 0.2      | 🟢      |
| Multigraphs      | 2%    | STRidge   | 10.0    | 12 | no  | 27/0  | 16       | 0.2379  | 0.1227 | 0.2      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | no  | 1/0   | 4        | 0.3372  | 0.1866 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9984  | 0.2228 | 0.2      | 🔴      |
| Multigraphs      | 2%    | STRidge   | 1.0     | 4  | yes | 22/0  | 14       | 1.0008  | 0.3485 | 0.3      | 🔴      |
| Multigraphs      | 2%    | MSTLS     | 0.25    | 12 | yes | 1/0   | 4        | 0.3329  | 0.1858 | 0.3      | 🔴      |
| Multigraphs      | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9984  | 0.2228 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 0.001   | 12 | no  | 37/0  | 20       | 0.8344  | 0.3033 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | no  | 5/0   | 5        | 0.9730  | 0.3915 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | no  | 0/5   | 1        | 0.9988  | 0.4485 | 0.2      | 🔴      |
| Multigraphs      | 5%    | STRidge   | 1.0     | 4  | yes | 30/0  | 17       | 1.0002  | 0.5244 | 0.3      | 🔴      |
| Multigraphs      | 5%    | MSTLS     | 0.75    | 8  | yes | 5/0   | 5        | 0.9617  | 0.3901 | 0.4      | 🔴      |
| Multigraphs      | 5%    | RFE       |         | 8  | yes | 0/5   | 1        | 0.9988  | 0.4485 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | no  | 0/0   | 3        | 0.0024  | 0.0112 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | no  | 0/0   | 3        | 0.0022  | 0.0112 | 0.2      | 🟢      |
| Mod. Multigraphs | Clean | STRidge   | 10.0    | 12 | yes | 1/0   | 4        | 0.0296  | 0.0100 | 0.2      | 🔴      |
| Mod. Multigraphs | Clean | MSTLS     | 0.5     | 12 | yes | 0/0   | 3        | 0.0278  | 0.0158 | 0.3      | 🟢      |
| Mod. Multigraphs | Clean | RFE       |         | 12 | yes | 0/0   | 3        | 0.0278  | 0.0158 | 0.2      | 🟢      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | no  | 27/0  | 16       | 0.2379  | 0.1227 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | no  | 1/0   | 4        | 0.3372  | 0.1866 | 0.4      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9984  | 0.2228 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | STRidge   | 10.0    | 12 | yes | 34/0  | 17       | 0.4032  | 0.1233 | 0.2      | 🔴      |
| Mod. Multigraphs | 2%    | MSTLS     | 0.25    | 12 | yes | 1/0   | 4        | 0.3329  | 0.1858 | 0.3      | 🔴      |
| Mod. Multigraphs | 2%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9984  | 0.2228 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 0.1     | 4  | no  | 31/0  | 18       | 1.0004  | 0.5244 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.75    | 12 | no  | 10/0  | 8        | 0.8007  | 0.3099 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | no  | 0/5   | 1        | 0.9988  | 0.3690 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | STRidge   | 1.0     | 12 | yes | 30/0  | 16       | 0.9072  | 0.3098 | 0.2      | 🔴      |
| Mod. Multigraphs | 5%    | MSTLS     | 0.25    | 12 | yes | 0/5   | 1        | 0.9988  | 0.3690 | 0.3      | 🔴      |
| Mod. Multigraphs | 5%    | RFE       |         | 12 | yes | 0/5   | 1        | 0.9988  | 0.3690 | 0.2      | 🔴      |


## Oracle (ground-truth residual)

The oracle plugs the true coefficients into ``evaluate_unified`` (same pipeline as all regression results).  Oracle E = 0 always (w ≡ true).  E_fit = ‖Θ·true − Ft‖/‖Ft‖ is the irreducible noise floor.

| Noise | Formulation   | Backend    | m  | IBP | E_data floor | cond(Θ)     |
|:------|:--------------|:-----------|:--:|:---:|-------------:|:------------|
| Clean | weak (no IBP) | spectral   | 4  | no  | 0.0171       | 15559622.6  |
| Clean | weak (no IBP) | spectral   | 8  | no  | 0.0142       | 21531316.4  |
| Clean | weak (no IBP) | spectral   | 12 | no  | 0.0114       | 34013918.5  |
| Clean | weak (no IBP) | polynomial | 4  | no  | 0.0227       | 15564349.3  |
| Clean | weak (no IBP) | polynomial | 8  | no  | 0.0145       | 21542667.8  |
| Clean | weak (no IBP) | polynomial | 12 | no  | 0.0114       | 34063590.5  |
| Clean | weak (IBP)    | spectral   | 4  | yes | 178.9467     | 59713062.9  |
| Clean | weak (IBP)    | spectral   | 8  | yes | 0.7600       | 22310889.7  |
| Clean | weak (IBP)    | spectral   | 12 | yes | 0.0326       | 34140137.1  |
| Clean | weak (IBP)    | polynomial | 4  | yes | 178.9467     | 59112570.3  |
| Clean | weak (IBP)    | polynomial | 8  | yes | 0.7599       | 22378946.1  |
| Clean | weak (IBP)    | polynomial | 12 | yes | 0.0326       | 34199127.9  |
| 2%    | weak (no IBP) | spectral   | 4  | no  | 18.4247      | 21474809.9  |
| 2%    | weak (no IBP) | spectral   | 8  | no  | 11.3180      | 33213146.7  |
| 2%    | weak (no IBP) | spectral   | 12 | no  | 11.4029      | 45377815.8  |
| 2%    | weak (no IBP) | polynomial | 4  | no  | 5.3848       | 18300991.0  |
| 2%    | weak (no IBP) | polynomial | 8  | no  | 0.9636       | 24727456.6  |
| 2%    | weak (no IBP) | polynomial | 12 | no  | 0.8842       | 32039332.3  |
| 2%    | weak (IBP)    | spectral   | 4  | yes | 186.5849     | 21421560.4  |
| 2%    | weak (IBP)    | spectral   | 8  | yes | 11.9736      | 33218011.5  |
| 2%    | weak (IBP)    | spectral   | 12 | yes | 11.4293      | 45383276.7  |
| 2%    | weak (IBP)    | polynomial | 4  | yes | 177.7529     | 18291220.6  |
| 2%    | weak (IBP)    | polynomial | 8  | yes | 1.6306       | 24702480.8  |
| 2%    | weak (IBP)    | polynomial | 12 | yes | 0.9103       | 32035475.0  |
| 5%    | weak (no IBP) | spectral   | 4  | no  | 75.2191      | 64231689.3  |
| 5%    | weak (no IBP) | spectral   | 8  | no  | 69.4094      | 93763895.2  |
| 5%    | weak (no IBP) | spectral   | 12 | no  | 70.6711      | 125136634.2 |
| 5%    | weak (no IBP) | polynomial | 4  | no  | 13.4584      | 18973257.7  |
| 5%    | weak (no IBP) | polynomial | 8  | no  | 5.4893       | 27941601.0  |
| 5%    | weak (no IBP) | polynomial | 12 | no  | 5.4493       | 36171050.9  |
| 5%    | weak (IBP)    | spectral   | 4  | yes | 226.1469     | 63892845.2  |
| 5%    | weak (IBP)    | spectral   | 8  | yes | 70.0482      | 93755748.7  |
| 5%    | weak (IBP)    | spectral   | 12 | yes | 70.6973      | 125136813.2 |
| 5%    | weak (IBP)    | polynomial | 4  | yes | 172.5340     | 18941693.6  |
| 5%    | weak (IBP)    | polynomial | 8  | yes | 6.1362       | 27917652.5  |
| 5%    | weak (IBP)    | polynomial | 12 | yes | 5.4753       | 36169050.1  |
---
| Source                                        | Wall time | CPU                                           | Threads |
|:----------------------------------------------|----------:|:----------------------------------------------|:-------:|
| `sweep_20260715_190818_ks_wchaos_trimmed.csv` | 45.7 min  | AMD Ryzen 7 PRO 7840U w/ Radeon 780M Graphics | 1       |
