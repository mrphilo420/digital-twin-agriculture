# Multivariate Lattice Digital Twin for Agricultural Soil Nutrient Dynamics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-≥4.0.0-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

A spatially explicit Digital Twin framework for modeling multivariate (N-P-K) soil nutrient dynamics under crop rotation systems. This implementation validates theoretical results on mass conservation, spatial autocorrelation, and variance scaling under α-mixing dependence structures.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Simulation Scenarios](#simulation-scenarios)
- [Results Summary](#results-summary)
- [Theoretical Proofs](#theoretical-proofs)
- [Repository Structure](#repository-structure)
- [Citation](#citation)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## 📖 Overview

This repository implements a **Multivariate Lattice Digital Twin** for agricultural soil nutrient management. The framework addresses three critical gaps in existing approaches:

1. **Spatial heterogeneity**: Conventional soil sampling aggregates to field-averaged values, missing fine-scale variation
2. **Multivariate interactions**: Single-nutrient analyses ignore N-P-K coupling and competition
3. **Spatial autocorrelation**: Classical statistical tests assume independence, inflating Type I error rates

The Digital Twin represents agricultural soil as a 4-D tensor (x, y, time, nutrient channel) with:
- Spatially variable buffering capacity (α) based on soil texture
- Crop-specific nutrient removal forces from empirical literature
- FFT-based Gaussian kernel smoothing for lateral transport
- Proven mass conservation and positivity preservation

**Manuscript Status**: Under review at *European Journal of Agronomy* (Elsevier)

---

## ✨ Key Features

### Methodological Contributions

| Feature | Description | Validation |
|---------|-------------|------------|
| **Mass Conservation** | Reflective-boundary Gaussian smoothing with renormalization | Theorem 1 |
| **Spatial Autocorrelation** | Moran's I via O(N) pure-summation algorithm | Theorem 2 |
| **Variance Scaling** | Stress variance ∝ buffering variance (short horizons) | Proposition 1 |
| **Positivity Preservation** | Multiplicative update guarantees S > 0 | Proposition 2 |
| **Dynamic Scaling Law** | R(T) increases monotonically with time | Lemma 2 (extended) |

### Agronomic Insights

- **Diverse rotations** (corn-soybean-wheat) preserve substantial phosphorus stress signal (38% of total variance)
- **Continuous corn** collapses P-stress to 7.7%, masking progressive depletion
- **Sandy zones** (low buffering) show 2× higher depletion than clay-rich areas
- **Spatial structure** remains stable (Moran's I = 0.998) across rotation scenarios

---

## 📦 Installation

### Prerequisites

- R version ≥ 4.0.0
- Git (for cloning the repository)

### Required R Packages

```r
# Core simulation packages
install.packages(c("fields", "spdep", "ggplot2", "patchwork", 
                   "viridis", "scales", "RColorBrewer", 
                   "reshape2", "dplyr", "tidyr", "cowplot"))

# Optional: for parallel computation
install.packages("doParallel")
```

### Clone Repository

```bash
git clone https://github.com/mrphilo420/digital-twin-agriculture.git
cd digital-twin-agriculture
```

---

## 🚀 Quick Start

### Run Default Simulation

```r
# Load main simulation script
source("src/run_simulation.R")

# This will:
# 1. Create buffering capacity map (100x100 grid)
# 2. Run baseline rotation (corn-soybean-wheat, 3 years)
# 3. Run continuous corn (3 years)
# 4. Generate summary statistics
# 5. Save results to results/ directory
```

### Generate Figures

```r
# Create journal-compliant figures
source("src/generate_figures.R")

# Output: figures/ directory with TIFF (600 DPI), PDF, and PNG files
```

### Reproduce Key Results

```r
# Load saved results
baseline <- readRDS("results/baseline_rotation.rds")
continuous <- readRDS("results/continuous_corn.rds")

# Print stress decomposition
cat("Baseline rotation stress weights:\n")
print(baseline$stress_weights)

cat("\nContinuous corn stress weights:\n")
print(continuous$stress_weights)

# Verify Moran's I
cat("\nMoran's I (Baseline):", baseline$moran_I, "\n")
cat("Moran's I (Continuous):", continuous$moran_I, "\n")
```

---

## 🔬 Simulation Scenarios

### Default Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Grid size | 100 × 100 | 10,000 spatial locations |
| Resolution | Δ = 0.01 | Normalized units (~10 m) |
| Time horizon | T = 3 years | Standard rotation cycle |
| Extended | T = 15 years | Long-run depletion analysis |
| Kernel σ | 1.2 cells | Gaussian smoothing bandwidth |
| Nutrients | N, P, K | Three-channel tensor |

### Rotation Scenarios

1. **Baseline**: Corn → Soybean → Wheat (3 years)
2. **Continuous Corn**: Corn → Corn → Corn (3 years)
3. **Null**: Zero forces (mass conservation test)

### Crop Force Vectors

| Crop | N | P | K | Source |
|------|---|---|---|--------|
| Corn | -0.60 | -0.20 | -0.20 | You et al. (2023) |
| Soybean | +0.20 | -0.10 | -0.10 | Salvagiotti et al. (2008) |
| Wheat | -0.20 | -0.40 | -0.10 | Liu et al. (2011) |

*Note: Positive N force for soybean represents biological fixation*

---

## 📊 Results Summary

### Stress Decomposition (3-Year)

| Scenario | Mean D | Max D | w_N | w_P | w_K | Moran's I |
|----------|--------|-------|-----|-----|-----|-----------|
| **Baseline** | 0.516 | 0.909 | 56.2% | 38.0% | 5.7% | 0.998 |
| **Continuous Corn** | 0.762 | 1.163 | 84.6% | 7.7% | 7.7% | 0.998 |

### Lemma 2 Variance Scaling

| Time Horizon | R = Var[D]/Var[α] | Interpretation |
|--------------|-------------------|----------------|
| T = 3 years | 0.769 | First-order approximation valid |
| T = 15 years | 1.231 | Non-linear amplification evident |

### Long-Run Depletion Trajectory (Baseline)

| Cycle | Year | Mean N (fraction) |
|-------|------|-------------------|
| 1 | 3 | 0.667 |
| 2 | 6 | 0.490 |
| 3 | 9 | 0.375 |
| 4 | 12 | 0.293 |
| 5 | 15 | 0.213 |

*No steady state observed; monotonic decline toward zero (agronomically correct)*

---

## 📐 Theoretical Proofs

The Appendix contains formal proofs for:

1. **Theorem 1**: Mass conservation of renormalized kernel smoothing
2. **Theorem 2**: Equivalence of summation and matrix forms of Moran's I
3. **Proposition 1**: First-order variance scaling (Var[D] ∝ Var[α])
4. **Proposition 2**: Positivity preservation (S_{t,c} > 0 ∀ t)
5. **Theorem 3**: Weak convergence of smoothed empirical process
6. **Theorem 4**: Asymptotic null distribution of spatial CvM statistic

See `manuscript/appendix_proofs.pdf` for complete mathematical derivations.

---

## 📁 Repository Structure

```
digital-twin-agriculture/
├── README.md                      # This file
├── LICENSE                        # MIT License
├── CITATION.cff                   # Citation metadata
│
├── src/                           # Source code
│   ├── run_simulation.R           # Main simulation driver
│   ├── create_buffering_map.R     # Alpha map generation
│   ├── gaussian_smooth.R          # Kernel smoothing functions
│   ├── calc_morans_I.R            # Spatial autocorrelation
│   └── generate_figures.R         # Visualization scripts
│
├── results/                       # Simulation outputs (generated)
│   ├── baseline_rotation.rds      # Baseline scenario results
│   ├── continuous_corn.rds        # Monoculture results
│   ├── alpha_map.rds              # Buffering capacity map
│   └── summary_statistics.csv     # Key metrics table
│
├── figures/                       # Journal-ready figures (generated)
│   ├── fig1_buffering_capacity.tiff
│   ├── fig2_stress_comparison.tiff
│   ├── fig3_stress_decomposition.tiff
│   ├── fig4_morans_I.tiff
│   └── graphical_abstract.tiff
│
├── manuscript/                    # Paper files
│   ├── main_text.tex              # LaTeX manuscript
│   ├── appendix_proofs.pdf        # Mathematical proofs
│   ├── references.bib             # Bibliography
│   └── response_to_reviewers.pdf  # Revision history
│
├── data/                          # External data (if applicable)
│   └── (none - all data synthetic)
│
└── tests/                         # Unit tests
    ├── test_mass_conservation.R
    ├── test_positivity.R
    └── test_morans_I.R
```

---

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{mandap2026digital,
  title = {Multivariate Lattice Deformation: A Spatially Explicit Digital Twin for Predicting Nutrient Dynamics in Heterogeneous Crop Rotations},
  author = {Mandap, Marco},
  journal = {European Journal of Agronomy},
  year = {2026},
  volume = {under review},
  doi = {10.xxxx/xxxxx},
  url = {https://github.com/mrphilo420/digital-twin-agriculture}
}
```

**Preprint**: Available on [SSRN/arXiv] (link pending)

---

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

**Summary**: You are free to use, modify, and distribute this code for academic and commercial purposes, provided you:
- Include the original copyright notice
- Include the license text
- Acknowledge the original authors

---

## 🙏 Acknowledgments

### Funding
This research received no specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

### Computational Resources
Simulations were performed on Bulacan State University using R version 4.4.3.

### Software Dependencies
- **fields**: Nychka et al. (2021) - Spatial analysis and FFT smoothing
- **spdep**: Bivand et al. (2023) - Spatial weights and Moran's I
- **ggplot2**: Wickham (2016) - Visualization framework
- **patchwork**: Pedersen (2023) - Multi-panel figure composition

### Development History

| Version | Date | Changes |
|---------|------|---------|
| v1.0 | Jan 2026 | Initial implementation (5-D tensor, depth layers) |
| v2.0 | Feb 2026 | Fixed MM kinetics for N fixation (biological inconsistency) |
| v3.0 | Mar 2026 | Corrected Lemma 2 direction (Var[α] not Var[1/α]) |
| v4.0 | Apr 2026 | Simplified to 4-D formulation; reflective boundaries |
| v4.1 | May 2026 | Added extended simulation (T=15); falsified steady-state conjecture |

**Development Log**: See `docs/development_log.md` for detailed revision history.

---

## 📧 Contact

**Corresponding Author**: [Marco Mandap, PhD]
**Email**: [marco.mandap@bulsu.edu.ph]  
**GitHub**: [@mrphilo420](https://github.com/mrphilo420)  

**Issues & Bug Reports**: Please use the [GitHub Issues](https://github.com/mrphilo420/digital-twin-agriculture/issues) tracker.

---

## 🔗 Related Resources

- **Elsevier LaTeX Template**: https://www.ctan.org/pkg/elsarticle
- **R Spatial Task View**: https://cran.r-project.org/web/views/Spatial.html
- **Digital Twin in Agriculture Review**: Purcell & Neubauer (2023)
- **Spatial Statistics Textbook**: Cressie (1993), Goovaerts (1997)

---

<div align="center">
  
**Last Updated**: May 2026  
**Status**: ✅ Active Development | 📝 Under Review
