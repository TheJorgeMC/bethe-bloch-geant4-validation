# Theoretical-Computational Validation of the Bethe-Bloch Equation for Protons in the Clinical Proton Therapy Range

Theoretical derivation, analytic implementation, and Monte Carlo validation of the Bethe-Bloch stopping power formula for protons, cross-checked against NIST reference data.

🚧 **Status:** active development (early stage).

## Overview

This project validates the Bethe-Bloch equation for protons through three independent approaches:

1. **Theory** — full derivation of the Bethe-Bloch formula (relativistic, with shell and density corrections), implemented analytically in both **Python and MATLAB** (cross-validated against each other).
2. **Monte Carlo simulation** — particle transport simulation in **Geant4**.
3. **Reference data** — stopping power tables from **NIST PSTAR**.

The three results are compared against each other across an energy range of **3–300 MeV**, matching the clinically relevant proton therapy energy range (Paganetti, 2012), for four absorber materials spanning a wide range of atomic number: **water, aluminum, copper, and lead**. Particular attention is given to the energy regions where the Bethe-Bloch approximation is expected to break down (low energies, below βγ≈0.1, i.e. T≲4.7 MeV for protons), which is directly relevant to proton therapy applications near the Bragg peak.

## Repository structure

```
theory/             Full theoretical derivation (LaTeX/Markdown)
python/
  analytic/         Python implementation of the Bethe-Bloch formula
  analysis/         Comparison, plotting, and error analysis (Python)
matlab/
  analytic/         MATLAB implementation of the Bethe-Bloch formula
  analysis/         Comparison, plotting, and error analysis (MATLAB)
geant4/             Geant4 simulation source code
nist_data/          NIST PSTAR reference tables (parsed to CSV)
figures/            Generated plots
report/             Final technical report (PDF)
```

## Methodology

| Source | Method |
|---|---|
| Theory | Analytic Bethe-Bloch formula (Python + MATLAB, cross-validated) |
| Simulation | Geant4 Monte Carlo, monoenergetic proton beams through thin absorber slabs |
| Reference | NIST PSTAR stopping power tables |

## Requirements

- **Geant4** (LTS release recommended)
- **Python 3.x** — numpy, scipy, pandas, matplotlib
- **MATLAB** — base installation, no specialized toolboxes required

Detailed setup and run instructions will be added as each module is completed.

## Motivation

Proton stopping power is central to proton therapy treatment planning, where the depth-dose profile (and the location of the Bragg peak) depends directly on how accurately energy loss can be predicted. The 3–300 MeV energy range studied here matches the clinical range used in proton therapy (Paganetti, 2012), and was chosen deliberately to span both the region where the standard Bethe-Bloch theory holds (βγ≳0.1) and the region near its lower validity limit, where shell and Barkas corrections become significant.

## References

- Paganetti, H. (2012). *Proton Therapy Physics*. CRC Press.
- Particle Data Group, "Passage of Particles Through Matter," *Prog. Theor. Exp. Phys.*

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Author

[TheJorgeMC](https://github.com/TheJorgeMC)
