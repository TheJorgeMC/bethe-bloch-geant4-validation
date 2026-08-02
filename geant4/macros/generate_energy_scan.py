#!/usr/bin/env python3
# ============================================================================
# generate_energy_scan.py — single source of truth for the energy sweep.
#
# From one energy grid (by default the NIST PSTAR preferred-number grid read
# from ../../nist_data/energy_grid.csv, for direct cross-validation against
# PSTAR) this script generates BOTH sweep mechanisms, so they can never drift
# apart:
#
#   1) run_scan/run_<Etag>MeV.mac — one archivable macro per curve point
#      (data <-> conditions traceability, design section 7.2).
#   2) energy_scan.mac (written next to this script) — a single macro with
#      the per-point commands INLINED (no /control/foreach: see note below).
#
# THIN-SLAB THICKNESS RULE (the key methodological point): the slab
# thickness is chosen PER ENERGY as
#
#     t(E) = clamp( FRAC_RANGE x R(E),  T_MIN,  T_MAX )
#
# with R(E) the proton CSDA range in water. A fixed thickness silently
# changes the nature of the measurement along the sweep: 5 mm is a thin
# slab above ~100 MeV but exceeds the proton range below ~20.5 MeV, where
# every primary stops inside and the total (~unrestricted) estimator
# (E_in - E_out)/track — the ONLY one comparable with the Bethe formula
# and with PSTAR — is undefined. Scaling t with E keeps the fractional
# energy loss ~FRAC_RANGE at every point, so simulation, theory and
# tabulated data can be validated over the SAME energy range with no
# systematic gaps. R(E) uses the Bragg-Kleeman rule R = alpha x E^p for
# protons in water, with alpha = 0.002777 cm/MeV^p and p = 1.723 under the
# CSDA approximation, taken from doi:10.48550/arXiv.2011.00285 — the
# documented basis for both the rule and its constants. Any residual
# inaccuracy of the fit only shifts the ~5% target loss slightly and is
# clamped by T_MAX; it plays no role in the physics results (the analysis
# divides by the actual per-event track length). Below ~0.3 MeV even T_MIN
# exceeds the range; those grid points remain stopping measurements and
# are excluded from the dE/dx validation band (the analysis flags them via
# n_exit = 0).
#
# FILE-NAME TAGS: output file names are now
# "dedx_<materialKey>_<Etag>MeV_<cutTag>" (e.g. "dedx_water_150MeV_cut0p01mm"),
# so analyze_dedx.C can auto-detect which material a run belongs to instead
# of assuming water. G4AnalysisManager treats everything after the last '.'
# as a file-type extension, so decimal energies like 0.001 are 'p'-encoded:
# 0.001 MeV -> "dedx_water_0p001MeV_cut0p01mm". The /gun/energy command
# keeps the numeric value. Same reason there is no /control/foreach: each
# point needs the numeric energy AND the 'p'-tag (plus now a per-point
# thickness), and a single foreach alias cannot carry them all — points are
# written inline. run_point.mac remains only as a deprecated stub for the
# CMake configure_file() rule.
#
# BREAKING CHANGE: file names before this version had no material tag
# (just "dedx_150MeV_cut0p01mm..."). Existing water CSV/mac files must
# either be renamed to insert "water_" after "dedx_" (or "run_water_" for
# the run_scan/ macros) or be re-simulated; analyze_dedx.C's material
# auto-detection cannot recognize the old, untagged names.
#
# Usage:
#   python3 generate_energy_scan.py            # NIST PSTAR grid
#   python3 generate_energy_scan.py --log 20   # 20 log points in [E_MIN, E_MAX]
#
# Then either:
#   for m in run_scan/run_*MeV.mac; do ./slab "$m"; done
# or:
#   ./slab macros/energy_scan.mac
# ============================================================================
import argparse
import math
import pathlib
import numpy as np
import pandas as pd
import hashlib

import bragg_kleeman_materials as bk

# --- Material -----------------------------------------------------------
# One of bk.MATERIALS: "water", "air", "aluminium", "copper", "lead".
# analyze_dedx.C auto-detects the material of a run from the tag embedded
# in the file name below (dedx_<MATERIAL_KEY>_...), so switching this and
# re-running the sweep is enough to retarget the whole analysis -- EXCEPT
# for "air", which analyze_dedx.C does not recognize (no sourced mean
# excitation energy I for the Bethe formula; air is only a bridging
# material for the Bragg-Kleeman range-scaling rule, not a dE/dx
# validation target). See referencias_y_decisiones.md.
MATERIAL_KEY = "water"
MATERIAL = bk.G4_MATERIAL_NAME[MATERIAL_KEY]  # Geant4 NIST material name

# Materials included in the single energy_scan.mac sweep.
# "air" is deliberately excluded: it's only the bridging material used by
# the Bragg-Kleeman derivation below, not a dE/dx validation target
# (analyze_dedx.C has no sourced mean excitation energy I for it -- see
# the MATERIAL_KEY comment above).
VALIDATION_MATERIALS = [k for k in bk.MATERIALS if k != "air"]

MATERIAL_DATA_DIR = pathlib.Path("material_data")

# --- Common sweep conditions (edit here, not in the macros) -----------------
SIZE_XY = "10 cm"
# Layers: 1 for the sweep (depth-dose segmentation is meaningless for
# micrometer-thin low-energy slabs; use run.mac with 50 layers for Bragg
# curve studies at a fixed energy).
N_LAYERS = 1
CUT = "0.01 mm"
CUT_TAG = "cut0p01mm"   # appears in the data file name
N_EVENTS = int(1e7)         # test statistics; raise for production

# --- Thin-slab thickness rule ------------------------------------------------
FRAC_RANGE = 0.05       # target fractional energy loss per slab (~5%)
T_MAX_MM = 5.0          # never thicker than the original design slab
T_MIN_MM = 0.001        # 1 um floor (buildability/step-size sanity)

# Bragg-Kleeman R = alpha * E^p for protons, CSDA approximation. Water's
# (alpha, p) come directly from doi:10.48550/arXiv.2011.00285; air,
# aluminium, copper and lead are DERIVED from water via the generalized
# Bragg-Kleeman range-scaling rule (Bragg & Kleeman, 1905), implemented in
# bragg_kleeman_materials.py -- see that module for the full derivation,
# the sourced material properties, and a self-check against the
# independently confirmed air reference values (rho = 1.2929e-3 g/cm3,
# sqrt(A_air) = 3.81). The fit grid below is only used to recover
# (alpha, p) for aluminium/copper/lead by log-log regression through the
# range curve computed with the scaling formula -- it does NOT need to
# match the actual sweep grid, since the underlying relation is an exact
# power law (any sufficiently wide, well-sampled grid recovers the same
# constants to floating-point precision).
bk.verify_air_reference_values()
_FIT_ENERGIES_MEV = np.geomspace(3.0, 1000.0, 200)
_MATERIAL_CONSTANTS = bk.derive_all_constants(_FIT_ENERGIES_MEV)
BK_A_CM, BK_P = _MATERIAL_CONSTANTS[MATERIAL_KEY]


def slab_thickness_mm(E_MeV: float) -> float:
    """Thin-slab rule: 5% of the Bragg-Kleeman range, clamped to [1 um, 5 mm]."""
    range_cm = BK_A_CM * (E_MeV ** BK_P)
    t_mm = FRAC_RANGE * range_cm * 10.0
    return min(max(t_mm, T_MIN_MM), T_MAX_MM)

def slab_thickness_mm_for(material_key: str, E_MeV: float) -> float:
    """Same thin-slab rule as slab_thickness_mm(), but for any material in
    _MATERIAL_CONSTANTS -- not just the globally active MATERIAL_KEY."""
    alpha_cm, p = _MATERIAL_CONSTANTS[material_key]
    range_cm = alpha_cm * (E_MeV ** p)
    t_mm = FRAC_RANGE * range_cm * 10.0
    return min(max(t_mm, T_MIN_MM), T_MAX_MM)

# Energy grid: preferred-number array, same as NIST PSTAR data, for direct
# cross-validation. Identical to ../../nist_data/energy_grid.csv.
NIST_GRID_CSV = pathlib.Path("../../nist_data/energy_grid.csv")
E_MIN, E_MAX = 3.0, 1000.0  # range for --log

# Output locations: per-point macros go to run_scan/; energy_scan.mac and
# the run_point.mac stub are rewritten in place, next to this script.
OUT_DIR = pathlib.Path("run_scan")
MACRO_DIR = pathlib.Path(__file__).resolve().parent


def etag(e: float) -> str:
    """Energy tag for file names: %g formatting with '.' -> 'p'."""
    return f"{e:g}".replace(".", "p")

def seeds_for(material_key: str, e: float) -> tuple[int, int]:
    """
    Deterministic seed pair for a sweep point, derived
    ONLY from (material, energy) -- not from the grid order or size --
    so that regenerating macros (even with a different/enlarged grid)
    reproduces the SAME seed for each point that already existed.

    Values re limited well below 2^{31-1} = (2,147,483,647) to
    avoid a documented Geant4 bug where the parser of
                /random/setSeeds
    converts to a 32-bit int instead of a long, silently truncating
    large values (see https://geant4-forum.web.cern.ch/t/random-setseeds-seems-to-not-work/3541)
    """
    digest = hashlib.sha256(f"{material_key}_{e:.6f}".encode()).digest()
    seed1 = int.from_bytes(digest[0:4], "big")
    seed1 = int.from_bytes(digest[0:4], "big") % 900_000_000 + 1
    seed2 = int.from_bytes(digest[4:8], "big") % 900_000_000 + 1
    return seed1, seed2

POINT_TEMPLATE = """# Macro generated by generate_energy_scan.py — do NOT edit by hand.
# Sweep point: {energy:g} MeV, cut {cut} ({cut_tag})
# Thin-slab rule: t = {thickness:.6g} mm (~{frac:.0%} of the estimated range)
/control/verbose 0
/run/verbose 1
/event/verbose 0
/tracking/verbose 0

/absorber/material {material}
/absorber/thickness {thickness:.6g} mm
/absorber/sizeXY {size_xy}
/absorber/numberOfLayers {n_layers}

/physics/absorberCut {cut}

/run/initialize

/gun/particle proton
/gun/energy {energy:g} MeV

/analysis/setFileName dedx_{material_key}_{tag}MeV_{cut_tag}

/random/setSeeds {seed1} {seed2}
/run/beamOn {n_events}
"""

ENERGY_SCAN_HEADER = """# ============================================================================
# energy_scan.mac — energy AND material sweep for the dE/dx vs E curve
# GENERATED by generate_energy_scan.py — do NOT edit by hand; edit the
# constants (and/or the NIST grid CSV) in that script and re-run it.
# Usage: ./slab macros/energy_scan.mac
#
# Sweeps ALL of {materials_description} in a single process. Each material
# block re-issues /absorber/material and /run/initialize: the production
# cut is stored as a fixed RANGE ({cut}), but Geant4 converts it to a
# material-specific ENERGY threshold when it rebuilds the physics tables
# in /run/initialize -- so re-initializing per material is not optional
# (design, section 5.3).
#
# Energy grid: {grid_description} ({n_points} points per material).
# Thickness follows the thin-slab rule t(E) = clamp({frac:.0%} x range(E),
# {tmin:g} mm, {tmax:g} mm) using EACH material's own Bragg-Kleeman
# constants (bragg_kleeman_materials.py).
# Output files land in material_data/<material>/, created by this script.
# ============================================================================
/control/verbose 1
/run/verbose 1
/event/verbose 0
/tracking/verbose 0
# --- Common conditions (material-independent) ----------------------------
/absorber/sizeXY {size_xy}
/absorber/numberOfLayers {n_layers}
/physics/absorberCut {cut}
"""

ENERGY_SCAN_MATERIAL_HEADER = """
# ============================================================================
# Material: {material_key} ({material})
# ============================================================================
/absorber/material {material}
/run/initialize
/gun/particle proton
"""

ENERGY_SCAN_POINT = """
# --- {energy:g} MeV  (t = {thickness:.6g} mm, seed {seed1} {seed2})
/absorber/thickness {thickness:.6g} mm
/gun/energy {energy:g} MeV
/analysis/setFileName material_data/{material_key}/dedx_{material_key}_{tag}MeV_{cut_tag}
/random/setSeeds {seed1} {seed2}
/run/beamOn {n_events}
"""

RUN_POINT_STUB = """# ============================================================================
# run_point.mac — DEPRECATED stub, kept only so CMake's configure_file()
# rule for this path keeps working. The sweep uses inlined points in
# energy_scan.mac (numeric energy + 'p'-encoded file tag + per-point
# thin-slab thickness cannot ride on a single /control/foreach alias).
# Safe to delete once the CMakeLists.txt entry for macros/run_point.mac is
# removed.
# ============================================================================
"""


def load_energies(log_points: int) -> tuple[list[float], str]:
    """Return the energy list (MeV) and a short description of its origin."""
    if log_points > 1:
        step = (math.log10(E_MAX) - math.log10(E_MIN)) / (log_points - 1)
        energies = [round(10 ** (math.log10(E_MIN) + i * step), 1)
                    for i in range(log_points)]
        return energies, f"logarithmic, [{E_MIN:g}, {E_MAX:g}] MeV"
    if not NIST_GRID_CSV.is_file():
        raise SystemExit(
            f"error: {NIST_GRID_CSV} not found — run from the directory where "
            "that relative path is valid, or use --log N")
    energies = (pd.read_csv(NIST_GRID_CSV, sep=",", header=0)
                .iloc[:, 0].dropna().astype(float).to_list())
    return energies, f"NIST PSTAR grid from {NIST_GRID_CSV}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--log", type=int, metavar="N", default=0,
        help="use N logarithmically spaced points between "
             f"{E_MIN:g} and {E_MAX:g} MeV instead of the NIST grid")
    args = parser.parse_args()
    energies, grid_description = load_energies(args.log)
    # Guard against tag collisions (two grid values mapping to one file name).
    tags = [etag(e) for e in energies]
    if len(set(tags)) != len(tags):
        dupes = sorted({t for t in tags if tags.count(t) > 1})
        raise SystemExit(f"error: duplicate energy tags in the grid: {dupes}")

    # --- 1) One archivable macro per (material, energy) sweep point ---------
    # Mirrors energy_scan.mac exactly (same materials, same grid, same
    # seeds_for() calls) so every point energy_scan.mac ever runs also has
    # a standalone macro you can re-run in isolation -- for an outlier you
    # want to reproduce alone, or to resume a sweep that crashed partway
    # through without re-running the materials already completed.
    OUT_DIR.mkdir(exist_ok=True)
    n_stopping = 0
    for mk in VALIDATION_MATERIALS:
        alpha_cm, p = _MATERIAL_CONSTANTS[mk]
        for e in energies:
            t = slab_thickness_mm_for(mk, e)
            seed1, seed2 = seeds_for(mk, e)
            # Estimated range shorter than the thickness floor: the proton
            # will stop inside even the thinnest buildable slab.
            if alpha_cm * (e ** p) * 10.0 < T_MIN_MM:
                n_stopping += 1
            macro = POINT_TEMPLATE.format(
                energy=e, tag=etag(e), thickness=t, frac=FRAC_RANGE,
                material=bk.G4_MATERIAL_NAME[mk], material_key=mk,
                size_xy=SIZE_XY, n_layers=N_LAYERS, cut=CUT, cut_tag=CUT_TAG,
                n_events=N_EVENTS, seed1=seed1, seed2=seed2)
            path = OUT_DIR / f"run_{mk}_{etag(e)}MeV.mac"
            path.write_text(macro)
    print(f"wrote {len(VALIDATION_MATERIALS) * len(energies)} macros in "
          f"{OUT_DIR}/ ({len(VALIDATION_MATERIALS)} materials x "
          f"{len(energies)} points)")
    if n_stopping:
        print(f"note: {n_stopping} (material, energy) points have range < "
              f"{T_MIN_MM} mm floor — they remain stopping measurements "
              "(flagged by the analysis via n_exit = 0)")

    # --- 2) energy_scan.mac: same materials/grid, points inlined ------------
    for mk in VALIDATION_MATERIALS:
        (MATERIAL_DATA_DIR / mk).mkdir(parents=True, exist_ok=True)

    scan = ENERGY_SCAN_HEADER.format(
        materials_description=", ".join(VALIDATION_MATERIALS),
        grid_description=grid_description, n_points=len(energies),
        frac=FRAC_RANGE, tmin=T_MIN_MM, tmax=T_MAX_MM,
        size_xy=SIZE_XY, n_layers=N_LAYERS, cut=CUT)
    for mk in VALIDATION_MATERIALS:
        scan += ENERGY_SCAN_MATERIAL_HEADER.format(
            material_key=mk, material=bk.G4_MATERIAL_NAME[mk])
        for e in energies:
            t = slab_thickness_mm_for(mk, e)
            seed1, seed2 = seeds_for(mk, e)
            scan += ENERGY_SCAN_POINT.format(
                energy=e, tag=etag(e), thickness=t,
                material_key=mk, cut_tag=CUT_TAG, n_events=N_EVENTS,
                seed1=seed1, seed2=seed2)
    scan_path = MACRO_DIR / "energy_scan.mac"
    scan_path.write_text(scan)
    print(f"wrote {scan_path} ({len(VALIDATION_MATERIALS)} materials x "
          f"{len(energies)} points = "
          f"{len(VALIDATION_MATERIALS) * len(energies)} runs, "
          f"{grid_description})")
    excluded = sorted(set(bk.MATERIALS) - set(VALIDATION_MATERIALS))
    print(f"  materials: {', '.join(VALIDATION_MATERIALS)}"
          + (f"  (excluded: {', '.join(excluded)})" if excluded else ""))
    print(f"  output files under {MATERIAL_DATA_DIR}/<material>/")

    # --- 3) Deprecated run_point.mac stub (CMake compatibility) -------------
    point_path = MACRO_DIR / "run_point.mac"
    point_path.write_text(RUN_POINT_STUB)
    print(f"wrote {point_path} (deprecated stub)")
    print("\nRun the sweep with either:")
    print(f'  for m in {OUT_DIR}/run_*MeV.mac; do ./slab "$m"; done   '
          f'# one process per point, resumable/parallelizable')
    print("  ./slab macros/energy_scan.mac"
          "   # everything in one process")
    print("\nNOTE: re-run cmake (or copy the regenerated macros) so the "
          "build dir picks up the new energy_scan.mac.")
    print("\nBragg-Kleeman constants for all materials (water is the sourced "
          "anchor; the rest are derived — see bragg_kleeman_materials.py):")
    for key, (alpha, p) in _MATERIAL_CONSTANTS.items():
        marker = "  [swept]" if key in VALIDATION_MATERIALS else "  [not swept]"
        print(f"  {key:<12} alpha = {alpha:.6e} cm/MeV^p   p = {p:.4f}{marker}")


if __name__ == "__main__":
    main()