# Agricultural Digital Twin: 5D Spatiotemporal Soil Nutrient Engine

[![R Status](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This repository contains the official R implementation of the **Multivariate Lattice Deformation** framework, a production‑ready Agricultural Digital Twin for simulating spatial soil nutrient dynamics (N, P, K) under heterogeneous crop rotations.

Unlike standard field‑averaged agronomic models, this physics‑inspired engine operates on a high‑resolution 5D tensor (Space $X, Y, Z$, Time $T$, and Nutrients $C$). It bridges geostatistics, plant physiology, and computational physics to predict localized soil degradation and economic yield risk.

---

## Core Features (The 5 Theoretical Pillars)

| Pillar | Description | Mathematical Foundation |
|--------|-------------|--------------------------|
| **1. Michaelis‑Menten Kinetics** | Replaces static crop removal rates with dynamic feed‑forward biological loops, guaranteeing non‑negativity via an Exponential Euler scheme. | Theorem 2 (Boundedness) |
| **2. 3D Volumetric Lattices** | Simulates depth‑partitioned root extraction and vertical mass‑conserving leaching cascades (e.g., topsoil P accumulation vs. subsoil NO₃ leaching). | Theorem 3 (Mass Conservation) |
| **3. Anisotropic Transport** | Utilizes highly optimized Fast Fourier Transform (FFT) convolutions to simulate directional tillage and topographic runoff without edge‑artifact mass loss. | Theorem 3 & ADR convergence |
| **4. Economic Translation** | Converts theoretical physical stress into localized USD yield penalties using the synergistic Mitscherlich‑Baule crop response surface. | Theorem 4 (Strict Concavity) |
| **5. Ensemble Kalman Filter (EnKF)** | A memory‑safe subspace‑inversion data assimilation module that allows the Digital Twin to self‑correct using noisy satellite (e.g., NDVI) or proximal soil sensor data. | Theorem 5 (Optimality) |

---

## Installation & Dependencies

The engine is built entirely in R. It relies on `imager` (C++ backend) for massive FFT speedups, allowing it to easily scale to high‑resolution farm grids (e.g., $1000 \times 1000$ pixels).

```r
install.packages(c("ggplot2", "reshape2", "MASS", "xtable", "viridis", "imager"))
```

For the Ensemble Kalman Filter demonstration, the `MASS` package is required (already included).

---

## Quickstart

Clone the repository and run the master script. It will:

1. Set up a $200 \times 200 \times 3$ heterogeneous soil buffering map.
2. Execute a baseline 3‑year corn–soybean–wheat rotation.
3. Execute a continuous corn (monoculture) scenario for comparison.
4. Compute depth‑weighted nutrient availability for a subsequent corn crop.
5. Translate the nutrient state into economic loss (USD/ha) using the Mitscherlich‑Baule function.
6. Demonstrate data assimilation by updating the topsoil nitrogen field with synthetic satellite observations via the EnKF.
7. Generate publication‑ready heatmaps and LaTeX summary tables.

```bash
git clone https://github.com/mrphilo420/digital-twin-agriculture.git
cd digital-twin-agriculture
Rscript multivar-lattice.R
```

### Example Outputs

| File | Description |
|------|-------------|
| `Figure1_3D_Volumetrics.png` | Heatmaps revealing the divergence between topsoil (0–20 cm) and subsoil (40–60 cm) nitrogen after three years. |
| `Figure2_Economic_Loss.png` | Spatial map of USD/ha yield penalties for a subsequent corn crop. |
| `Figure_Master_DigitalTwin.png` | Combined figure (volumetrics + economics) for publication. |

Console output includes summary statistics (mean stress, max stress, Moran’s I, economic loss) and a LaTeX‑formatted table ready for inclusion in a manuscript.
---

## ⚠️ Performance Notes

- The current settings (`X = 200`, `Y = 200`, `Z = 3`) run comfortably on a typical laptop (approx. 30 seconds per scenario).
- For very large grids (> $500 \times 500$), the Moran’s I calculation in Section 8 becomes memory‑intensive. If you encounter memory issues, comment out that block or compute Moran’s I on a random subset of 1000 pixels.
- The EnKF implementation is already memory‑safe and works for grids of any size.

---

## 📄 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to open a pull request or an issue on GitHub.

---

## 🙏 Acknowledgments

This work builds on foundational literature in spatial statistics, numerical analysis, and agronomy, including the seminal work of Evensen (2003) on the Ensemble Kalman Filter, Ripley (1981) on spatial statistics, and Paris (1974) on the Mitscherlich‑Baule response surface. The FFT convolution routines are powered by the excellent `imager` package by Simon Barthelme.

---

**Happy modeling!** 🌱
