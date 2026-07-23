# slab — Bethe-Bloch validation with protons in Geant4

Geant4 (11.x) application that simulates monoenergetic protons crossing a
homogeneous absorber slab, to compare the simulated stopping power against
the analytical Bethe-Bloch equation. Skeleton based
on the `TestEm0` style.

## Build

```bash
mkdir build && cd build
cmake -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4 ..
make -j$(nproc)
```

## Run

```bash
./slab                          # interactive (runs macros/vis.mac)
./slab macros/run.mac           # batch, EM option4 (default)
./slab macros/run.mac 3         # batch, EM option3 (sensitivity check)
./slab macros/energy_scan.mac   # full energy sweep
```

`macros/energy_scan.mac` and the per-point macros under `run_scan/` are both
**generated** by `macros/generate_energy_scan.py` from a single energy grid
(by default the NIST PSTAR preferred-number grid, for direct cross-validation
against PSTAR), so the two sweep mechanisms below can never drift apart:

```bash
python3 macros/generate_energy_scan.py
for m in run_scan/run_*MeV.mac; do ./slab "$m"; done
# or, equivalently:
./slab macros/energy_scan.mac
```

Re-run the generator (and re-run/copy into the build dir) whenever the energy
grid or the sweep constants change — see "Thin-slab thickness rule" below for
why the slab thickness is not fixed across the sweep.

## Macro parameters

| Command | Default | Meaning |
|---|---|---|
| `/absorber/material <NIST>` | `G4_WATER` | slab material |
| `/absorber/thickness <v> <u>` | 5 mm | thickness along z (the energy sweep overrides this per point — see "Thin-slab thickness rule" below) |
| `/absorber/sizeXY <v> <u>` | 10 cm | side of the transverse cross section |
| `/absorber/numberOfLayers <N>` | 1 | N>1 segments the slab (Bragg curve) |
| `/physics/absorberCut <v> <u>` | 1 mm | production cut in `AbsorberRegion` |
| `/gun/energy <v> <u>` | 150 MeV | proton kinetic energy |
| `/analysis/setFileName <n>` | `slab_out` | output file (after `/run/initialize`) |

## Output

**Console** (end of run): mean and rms (straggling) of the energy deposit,
per-channel balance, simulated "local" dE/dx (~restricted) and "total" dE/dx
(E_inc − E_exit), and the `G4EmCalculator` reference values: restricted
dE/dx (with the active cut), **unrestricted** dE/dx (the one comparable with
textbook Bethe-Bloch) and CSDA range.

**Ntuple `slab`** (one row per event, MeV / mm):

| Column | Meaning |
|---|---|
| `E_incident` | initial kinetic energy of the primary |
| `Edep_primary` | primary's deposit in the slab (≈ restricted dE/dx × path) |
| `Edep_secondary` | deposit of delta rays and their progeny inside the slab |
| `E_escaped_secondary` | kinetic energy carried out of the slab by secondaries |
| `E_exit_primary` | primary energy on exit (−1 if it never exited) |
| `track_length_primary` | actual primary path inside the slab |

Key design check (section 5.3): per event,
`E_incident ≈ Edep_primary + Edep_secondary + E_escaped_secondary +
E_exit_primary`.

**Histogram `DepthEdep`**: energy deposited vs depth (Bragg curve when
`numberOfLayers > 1`), one bin per layer.

### File format and multi-threading

The default output is **CSV**. Geant4 has no native ntuple merging for CSV
in multi-threaded mode (each worker writes `*_nt_slab_t<i>.csv`), so at the
end of every run the application automatically merges all per-thread ntuple
files into a single `<name>_nt_slab.csv` (header kept once, data rows
concatenated — row order across threads is irrelevant since events are
independent) and deletes the per-thread files. Histograms are always merged
natively, so `<name>_h1_DepthEdep.csv` is a single file as well.

Reading the merged ntuple with pandas:

```python
import pandas as pd
cols = ["E_incident", "Edep_primary", "Edep_secondary",
        "E_escaped_secondary", "E_exit_primary", "track_length_primary"]
df = pd.read_csv("dedx_water_150MeV_cut1mm_nt_slab.csv",
                 comment="#", names=cols)
```

For ROOT output instead, use a `.root` extension in `/analysis/setFileName`
(e.g. `dedx_150MeV.root`): a single natively merged file, readable with
`uproot`.

## Physics notes

- The primary's continuous dE/dx is the **restricted** one, with `T_cut`
  set by `/physics/absorberCut` (translated to energy material by material);
  it depends on that choice and is not a physical constant. The
  **unrestricted** value is recovered by adding the energy of the explicit
  deltas (ntuple columns) or by lowering the cut until convergence
  (e.g. 0.01 mm).
- No hadronic physics (design decision): there are no inelastic nuclear
  reactions; "lost" primaries can only have stopped inside the slab.
- Always document the Geant4 version, physics list (option3/4) and
  production cut alongside every result.

## References and methodology notes

This section records the sourced constants and methodological decisions
behind the energy sweep and its analysis (`macros/generate_energy_scan.py`,
`analysis/analyze_dedx.C`), so results can be reproduced and re-derived
without digging through commit history. It mirrors the project's living
decisions log (`referencias_y_decisiones.md`).

### Thin-slab thickness rule (Bragg-Kleeman)

A **fixed** slab thickness cannot serve the whole energy sweep: 5 mm is thin
above ~100 MeV but exceeds the proton range below ~20.5 MeV, where every
primary stops inside the slab and the "total" (~unrestricted) estimator
`(E_incident - E_exit_primary) / track_length_primary` — the only one
comparable with the Bethe formula and with PSTAR — becomes undefined. That
gap is a systematic validation error, not an expected data-quality
condition: a simulation-vs-theory-vs-PSTAR comparison is only valid if all
three cover the *same* energy range.

`generate_energy_scan.py` fixes this by setting the thickness **per energy
point**:

```
t(E) = clamp( FRAC_RANGE * R(E),  T_MIN,  T_MAX )
R(E) = alpha * E^p            (Bragg-Kleeman rule, CSDA approximation)
```

with `FRAC_RANGE = 0.05` (target ~5% fractional energy loss per slab),
`T_MIN = 1 um`, `T_MAX = 5 mm`, and

- `alpha = 0.002777 cm/MeV^p`
- `p = 1.723`

for protons in water under the CSDA approximation, taken from
**doi:10.48550/arXiv.2011.00285**. Scaling the thickness with energy keeps
the total estimator defined across the whole validation band (`E >= 3 MeV`,
the same cutoff used in `analyze_dedx.C`); only the ~35 lowest-energy grid
points below ~0.3 MeV still fall under the 1 um floor and remain legitimate
stopping measurements, excluded from the dE/dx comparison and flagged by the
analysis via `n_exit = 0`.

### Analytic reference: pure relativistic Bethe, no corrections

`analyze_dedx.C` compares the simulation against the **relativistic Bethe
formula with no shell, Barkas, Bloch, or density-effect corrections** —
matching the project's `analytic_solution.ipynb` reference notebook exactly,
not the PDG exact-Tmax/Sternheimer form:

```
S(T) = K * z^2 * (Z/A) * (1/beta^2) * [ ln(2 * me*c^2 * beta^2*gamma^2 / I) - beta^2 ]
```

using the heavy-projectile approximation `Tmax ~= 2*me*c^2*beta^2*gamma^2`
absorbed into the single logarithm. Constants: CODATA 2022
(`me = 0.51099895069 MeV`, `Mp = 938.27208943 MeV`), `K = 0.307075
MeV*cm^2/mol`, water `Z/A = 0.555087` via the Bragg additivity rule
(H: w=0.111894, A=1.0080; O: w=0.888106, A=15.999), `I = 78 eV` (ICRU 90).
Any introduced-physics corrections from the Geant4 physics list itself
(option3/option4) are expected to make the simulation deviate slightly from
this uncorrected analytic curve; that deviation is part of what is being
validated, not a bug.

### Restricted vs. unrestricted dE/dx — which to compare

- **Restricted** (`Edep_primary / track_length_primary`): the primary's
  continuous energy loss below the production cut `T_cut`. Depends on the
  cut chosen with `/physics/absorberCut`, so it is **not** directly
  comparable with the Bethe formula or PSTAR.
- **Total / ~unrestricted** (`(E_incident - E_exit_primary) /
  track_length_primary`): requires the primary to exit the slab (guaranteed
  across the validation band by the thin-slab rule above). This is the
  quantity to compare against `bethe_no_corr_MeV_cm` and against PSTAR.

`dedx_summary.csv` reports both, plus `rel_err_bethe_vs_total_pct`, the
relative error of `bethe_no_corr_MeV_cm` against `S_total_MeV_cm2_g` (the
total/unrestricted estimator) — never against the restricted one.

### Straggling fit: Landau (x) Gauss convolution ("langaus")

The per-event energy-deposit histograms (`AnalyzeRun` in
`analyze_dedx.C`) are fit with the Landau (x) Gauss convolution
(`langaufun`, following ROOT's official `langaus.C`), not a plain Gaussian
or plain Landau. Neither pure shape matches the true (Vavilov) straggling
distribution well enough at high statistics: with tens of thousands of
events, per-bin statistical errors shrink and any systematic shape mismatch
dominates chi2/ndf, which grows roughly linearly with the event count. The
langaus convolution is the practical, statistically justified stand-in for
Vavilov across both thin (Landau-like) and thick (Gaussian-like) regimes;
goodness of fit is reported via `TMath::Prob(chi2, ndf)` rather than
chi2/ndf alone, since the latter is not comparable across different
statistics.

### File-naming convention (`p`-tag)

`G4AnalysisManager` treats everything after the last `.` in a file name as a
file-type extension, so decimal energies cannot appear literally in output
file names. Both `generate_energy_scan.py` (`etag()`) and `analyze_dedx.C`
(`EnergyTag()`) apply the same encoding — replace `.` with `p` — e.g.
`0.001 MeV -> dedx_0p001MeV_cut0p01mm`. `/gun/energy` still receives the
numeric value; only the output file name is tagged. This is also why the
sweep points are inlined in `energy_scan.mac` instead of using
`/control/foreach`: a single alias cannot carry the numeric energy, the
`p`-tag, and the per-point thin-slab thickness at once.