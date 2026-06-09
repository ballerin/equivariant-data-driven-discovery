# Experimental results — equation summaries

This file summarizes the main results across all the experiments. The details for each equation are in the linked files below.

| Equation | File | Status |
|---|---|---|
| Allen–Cahn | [allen_cahn.md](allen_cahn.md) | ✅ Exact — easy equation |
| Swift–Hohenberg | [swift_hohenberg.md](swift_hohenberg.md) |  |
| Cahn–Hilliard | [cahn_hilliard.md](cahn_hilliard.md) |  |
| Porous Medium | [porous_medium.md](porous_medium.md) |  |
| Kuramoto–Sivashinsky | [kuramoto_sivashinsky.md](kuramoto_sivashinsky.md) |  |
 

## Clean data — equivariant vs standard library (E_fit)

The E column for the equivariant library measures distance to ground truth;
for the standard library it measures the **irreducible fit error** — how
well the recovered PDE reproduces $f_t$, regardless of interpretability.

| Equation | Equiv E | Equiv terms | Std E_fit | Std terms | Verdict |
|---|---|---|---|---|---|
| Allen–Cahn | 0.005 | 3 | 0.003 | 58 | Fits data but 19× more terms |
| Swift–Hohenberg | 0.027 | 4 | 0.008 | 43 | Fits data, Δ²f splits across 3 cols |
| Cahn–Hilliard | 0.049 | 4 | **0.693** | 60 | 69% variance unexplained |
| KS cellular | 0.003 | 3 | **0.347** | 60 | 35% variance unexplained |
| KS weak chaotic | 0.003 | 3 | **0.348** | 60 | 35% variance unexplained |
| KS chaotic | 0.005 | 3 | **0.264** | 60 | 26% variance unexplained |
| Porous Medium | 0.001 | 2 | **0.109** | 59 | 11% variance unexplained |

The standard library fails in two distinct ways:

1. **Representable but not sparse** (AC, SH): fits $f_t$ well (E_fit ≈ 0.003–0.008)
   but needs 10–19× more columns because isotropic operators split across
   multiple partial-derivative columns.

2. **Structurally absent** (CH, KS, PME): cannot fit the data because
   $\|\nabla f\|^2 = f_x^2 + f_y^2$ is a product of two derivatives, and the
   standard library only has $f^p \cdot \partial^a f$ (single-derivative) columns.

## Noisy data

| Equation | 0% noise | 2% noise | 5% noise | Key method |
|---|---|---|---|---|
| Allen–Cahn | ✅ E=0.005 | ✅ E=0.005 | ✅ E=0.003 | WSINDy (all levels) |
| Swift–Hohenberg | ✅ E=0.027 | ✅ E=0.059 | ✅ E=0.33 | WSINDy IBP (noisy) |
| Cahn–Hilliard | ✅ E=0.049 | ❌ E=1.0 | ❌ E=1.0 | — |
| KS (all 3) | ✅ E≈0.003 | ❌ E=1.0 | ❌ E=1.0 | — |
| Porous Medium | ✅ E=0.001 | ✅ E=0.003 | ✅ E=0.14 | WSINDy (all levels) |

## Summary by equation

| Equation | File | Status |
|---|---|---|
| Allen–Cahn | [allen_cahn.md](allen_cahn.md) | ✅ Exact at all noise levels |
| Swift–Hohenberg | [swift_hohenberg.md](swift_hohenberg.md) | ✅ Exact at all noise levels (IBP at 5%) |
| Cahn–Hilliard | [cahn_hilliard.md](cahn_hilliard.md) | ✅ Clean only; noisy fails completely |
| Porous Medium | [porous_medium.md](porous_medium.md) | ✅ Exact at all noise levels |
| Kuramoto–Sivashinsky | [kuramoto_sivashinsky.md](kuramoto_sivashinsky.md) | ✅ Clean only; noisy fails completely |
