# Non-Parametric Density Estimation

 **Course project** — *Non-Parametric Inference and Sequential Models*  
 **Indian Statistical Institute**  


## Authors
| Name |
|---|
| Sourav Chakraborty | 
| Saswata Seal |
| Piwal Abhishek Satish |

---

## Overview

A comprehensive, simulation-driven study of non-parametric density estimation in R.
The project systematically explores the **Histogram Estimator**, the theoretical
foundations of the **Glivenko–Cantelli Theorem**, and **Kernel Density Estimation (KDE)**,
validated empirically across — 9 standard univariate distributions
and 9 two- or three-component mixture models.

---

## Repository Structure

```
NonParametric-Density-Estimation/
│
├── 00_setup.R                      # Packages, shared utilities, mixture objects
├── 01_histogram_origin.R           # Effect of variable origin X0 on histograms
├── 02_histogram_bandwidth.R        # Effect of variable bandwidth h on histograms
├── 03_histogram_samplesize.R       # Asymptotic behaviour as n grows
├── 04_histogram_variablebins.R     # Variable bin-width histograms
├── 05_glivenko_cantelli.R          # Empirical verification of the GC theorem
├── 06_kernel_density_estimation.R  # Full KDE analysis
└── README.md
```

---

## How to Run

Each script starts with `source("00_setup.R")`, which installs any missing
packages and creates all shared objects. Run scripts in order or individually:

```r
source("00_setup.R")
source("01_histogram_origin.R")
source("02_histogram_bandwidth.R")
source("03_histogram_samplesize.R")
source("04_histogram_variablebins.R")
source("05_glivenko_cantelli.R")
source("06_kernel_density_estimation.R")
```

### Dependencies

| Package
|---|---|
| `ggplot2` |
| `patchwork` | 
| `smoothmest` |

All packages are auto-installed by `00_setup.R` if not already present.

---

## Distributions Covered

### Standard Univariate
Normal · Exponential · Double Exponential (Laplace) · Uniform · Beta · Gamma · Weibull · Cauchy · Logistic

### Mixture Models
| # | Mixture |
|---|---|
| 1 | 0.5·N(−2,1) + 0.5·N(2,1) — Symmetric Normal |
| 2 | 0.3·N(0,1) + 0.7·N(4,0.5) — Asymmetric Normal |
| 3 | 0.3·N(−3,0.5) + 0.4·N(0,1) + 0.3·N(3,0.5) — Trinormal |
| 4 | 0.6·N(0,1) + 0.4·Exp(1) — Normal–Exponential |
| 5 | 0.5·N(5,1) + 0.5·Gamma(2,1) — Normal–Gamma |
| 6 | 0.7·Beta(2,5) + 0.3·U(0,1) — Beta–Uniform |
| 7 | 0.8·Logistic(0,1) + 0.2·Cauchy(0,1) — Logistic–Cauchy |
| 8 | 0.5·Weibull(2,1) + 0.5·Gamma(3,2) — Weibull–Gamma |
| 9 | 0.5·N(0,1) + 0.5·Cauchy(0,1) — Normal–Cauchy |

---

## Script Summaries

### `01_histogram_origin.R`
Varies the histogram anchor point X0 ∈ {0, 0.2, 0.4, 1, −0.15, −0.35, −0.45, −2}
with fixed h = 0.5 and n = 50.  Demonstrates that integer-multiple shifts of h produce
**identical** histograms while fractional shifts cause **peak splitting** that can
create false bimodality — worst for mixture and heavy-tailed distributions.

### `02_histogram_bandwidth.R`
Varies h ∈ {0.05, 0.19, 0.34, 0.48, 0.62, 0.76, 0.91, 1.00} with fixed X0 = 0 and n = 50.
Empirically shows the **bias–variance trade-off**: small h gives spiky/undersmoothed
histograms (high variance), large h gives flat/oversmoothed histograms (high bias).
For Cauchy, no single global h works for the sparse tails.

### `03_histogram_samplesize.R`
Varies n ∈ {10, 20, 30, 70, 100, 250, 500, 1000} with fixed h = 0.3 and X0 = 0.
Cauchy(0, 1) and mixtures with Cauchy components require large n and careful bandwidth
choice.

### `04_histogram_variablebins.R`
Compares three binning schemes for Normal, Cauchy, and Exponential:
- Structured adaptive grid (narrow centre, wide tails) — best for symmetric distributions
- Wide random grid from U(−10,10) — chaotic representation
- Positive random grid from U(0,5) — domain misspecification for symmetric distributions;
  incidentally good for Exponential whose support is [0, inf)

### `05_glivenko_cantelli.R`
Empirically tracks V_n = sup_x |f̂_n(x) − f(x)| for n up to 50,000 and demonstrates:
- **Convergence confirmed** when all theorem conditions hold (h(n) converges , uniformly
  continuous f, bounded-variation kernel)
- **Failure: uniform continuity** — Exp(1) and Uniform distributions cause V_n to saturate
- **Failure: bandwidth** — h_n = n^{−1/2} (diverging series) causes V_n to grow
- **Failure: kernel variation** — Oscillating Sine and Damped Hyper-Oscillating kernels
  cause V_n to saturate at a non-zero value

### `06_kernel_density_estimation.R`
Full KDE analysis with rug plots showing individual kernel bumps:
- All 9 standard and all 9 mixture distributions
- Bandwidth effect (h = 0.1 vs h = 0.7)
- Kernel comparison: Gaussian, Epanechnikov, Rectangular, Triangular — confirming
  all standard symmetric kernels produce nearly identical estimates; **bandwidth
  dominates over kernel choice**
- Beta kernel for **boundary-corrected** KDE on (0, 1) support

---

## Key Findings

| Finding | Script |
|---|---|
| Integer-multiple X0 shifts yield identical histograms; fractional shifts cause peak splitting | 01 |
| Bias–variance trade-off is controlled entirely by bandwidth h | 02 |
| Variance is inversely proportional to n; fixed-h bias persists even at n = 1000 | 03 |
| Adaptive bins only help when data-driven; random grids can worsen estimates | 04 |
| Uniform continuity of f is non-negotiable for V_n -> 0 | 05 |
| a convergent series h(n) achieves GC convergence; divergent series fails | 05 |
| Kernels of unbounded variation prevent uniform convergence | 05 |
| KDE is robust to kernel shape; bandwidth is the critical parameter | 06 |
| Beta kernels eliminate boundary bias for bounded-support distributions | 06 |

---

## Reference

A snapshot of the findings can be found in the : `presentation.pdf`

The full theoretical derivations and simulation results are documented in the
accompanying project report: `report.pdf`
