# 🌾 Agricultural Digital Twin: 5D Spatiotemporal Soil Nutrient Engine

[![R Status](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This repository contains the official R implementation of the **Multivariate Lattice Deformation** framework, a production‑ready Agricultural Digital Twin for simulating spatial soil nutrient dynamics (N, P, K) under heterogeneous crop rotations.

Unlike standard field‑averaged agronomic models, this physics‑inspired engine operates on a high‑resolution 5D tensor (Space $X, Y, Z$, Time $T$, and Nutrients $C$). It bridges geostatistics, plant physiology, and computational physics to predict localized soil degradation and economic yield risk.

---

## ✨ Core Features (The 5 Theoretical Pillars)

| Pillar | Description | Mathematical Foundation |
|--------|-------------|--------------------------|
| **1. Michaelis‑Menten Kinetics** | Replaces static crop removal rates with dynamic feed‑forward biological loops, guaranteeing non‑negativity via an Exponential Euler scheme. | Theorem 2 (Boundedness) |
| **2. 3D Volumetric Lattices** | Simulates depth‑partitioned root extraction and vertical mass‑conserving leaching cascades (e.g., topsoil P accumulation vs. subsoil NO₃ leaching). | Theorem 3 (Mass Conservation) |
| **3. Anisotropic Transport** | Utilizes highly optimized Fast Fourier Transform (FFT) convolutions to simulate directional tillage and topographic runoff without edge‑artifact mass loss. | Theorem 3 & ADR convergence |
| **4. Economic Translation** | Converts theoretical physical stress into localized USD yield penalties using the synergistic Mitscherlich‑Baule crop response surface. | Theorem 4 (Strict Concavity) |
| **5. Ensemble Kalman Filter (EnKF)** | A memory‑safe subspace‑inversion data assimilation module that allows the Digital Twin to self‑correct using noisy satellite (e.g., NDVI) or proximal soil sensor data. | Theorem 5 (Optimality) |

---

## 🚀 Installation & Dependencies

The engine is built entirely in R. It relies on `imager` (C++ backend) for massive FFT speedups, allowing it to easily scale to high‑resolution farm grids (e.g., $1000 \times 1000$ pixels).

```r
install.packages(c("ggplot2", "reshape2", "MASS", "xtable", "viridis", "imager"))
