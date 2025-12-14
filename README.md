# Energy-Efficient Sub-Connected Active RIS (MATLAB)

This repository presents a MATLAB-based reproduction and evaluation of an **energy-efficient sub-connected active Reconfigurable Intelligent Surface (RIS)** architecture for next-generation wireless communications.

The project reproduces and verifies the main conclusions of the paper:

> **K. Liu, Z. Zhang, L. Dai**,  
> *Active Reconfigurable Intelligent Surface: Fully-Connected or Sub-Connected?*  
> IEEE Communications Letters, 2022.

---

## Overview

Reconfigurable Intelligent Surfaces (RIS) enable programmable wireless propagation and are regarded as a key technology for **6G green communications**.  

However:
- **Passive RIS** suffers from the *multiplicative fading effect*
- **Fully-connected active RIS** provides strong performance but incurs **high power consumption and hardware complexity**

This project focuses on a **sub-connected active RIS architecture**, where:
- Multiple RIS elements share **one power amplifier**
- Each element maintains **independent phase control**

This design significantly reduces energy consumption while preserving most of the beamforming gain.

---

## Key Contributions

- Complete **MATLAB implementation** of sub-connected active RIS
- Energy efficiency (EE) maximization using **fractional programming**
- Joint optimization of:
  - BS beamforming
  - RIS phase shifts
  - RIS amplification coefficients
- Fair comparison among:
  - Passive RIS
  - Fully-connected active RIS
  - Sub-connected active RIS

---

## Methodology

- Downlink multi-user MISO system with an active RIS
- EE formulated as a **fractional optimization problem**
- Solved via:
  - **Dinkelbach transformation**
  - **Alternating optimization**
- Convex subproblems solved using **CVX**

The algorithm iterates until convergence and guarantees a stable solution.

---

## Simulation Setup

- Platform: **MATLAB + CVX**
- Typical parameters:
  - BS antennas: M = 6
  - RIS elements: N = 256
  - Users: K = 4
- Rician fading channels (3GPP-inspired)
- Monte Carlo simulations (500 realizations)
- Identical power constraints across all schemes

---

## Results Summary

Key reproduced results:

- **Sub-connected active RIS achieves ~21% higher energy efficiency** than fully-connected active RIS
- **Only ~11% spectral efficiency loss**
- Consistently outperforms:
  - Passive RIS
  - No-RIS baseline
- Moderate sub-array size (e.g., **T = 16**) offers the best EE–SE trade-off

These results closely match the findings reported in the original paper.

