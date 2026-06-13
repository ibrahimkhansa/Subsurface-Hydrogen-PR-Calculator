## Subsurface hydrogen Peng-Robinson Calculator

### Project Overview
This repository provides a semi-analytical thermodynamic screening tool for **Underground Hydrogen Storage (UHS)**. When storing hydrogen in depleted gas fields or saline aquifers, assuming ideal gas behavior results in massive capacity estimation errors. Hydrogen features a unique molecular structure that strongly resists compression under high subsurface pressures, requiring a rigorous cubic Equation of State (EoS) to model correctly.

### Features
* **Multi-Component Solver:** Evaluates pure components ($H_2$, $CH_4$, $CO_2$) and dynamic blended streams.
* **Thermodynamic Precision:** Solves the cubic Peng-Robinson EoS from first principles to calculate the gas deviation factor ($Z$).
* **Subsurface Volumetrics:** Computes the Gas Formation Volume Factor ($B_g$) to translate downhole pore volumes into surface volumes (Standard Cubic Feet).
* **Advanced Mixing:** Employs quadratic van der Waals mixing rules to characterize gas-plume mixing zones.

### Core Engineering Insights
When running `generate_plots.py`, the calculator demonstrates a critical subsurface trade-off:
1. **Methane ($CH_4$)** compresses easily due to intermolecular attractions, dropping significantly below $Z = 1.0$.
2. **Hydrogen ($H_2$)** exhibits an aggressive resistance to compression under deep reservoir pressures, climbing well above $Z = 1.0$. 
3. **Mixtures** trace a unique fluid profile, proving that multi-component thermodynamic tracking is vital for accurate asset planning and site screening.
