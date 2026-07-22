# slab — Bethe-Bloch validation with protons in Geant4

Geant4 (11.x) application that simulates monoenergetic protons crossing a
homogeneous absorber slab, to compare the simulated stopping power against
the analytical Bethe-Bloch equation (Phase 5 of the project). Skeleton based
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

Alternative to the `foreach` sweep (one archivable macro per point):

```bash
python3 macros/generate_energy_scan.py
for m in run_scan/run_*MeV.mac; do ./slab "$m"; done
```

## Macro parameters

| Command | Default | Meaning |
|---|---|---|
| `/absorber/material <NIST>` | `G4_WATER` | slab material |
| `/absorber/thickness <v> <u>` | 5 mm | thickness along z |
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
