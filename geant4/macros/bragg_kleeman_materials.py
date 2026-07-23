#!/usr/bin/env python3
# ============================================================================
# bragg_kleeman_materials.py — generalized Bragg-Kleeman range rule for
# protons in water, air, aluminium, copper and lead.
#
# TWO DIFFERENT "Bragg" RULES are used side by side in this project and
# must not be confused:
#
#   1) The Bragg ADDITIVITY rule for Z/A (ICRU 37/49, Sec. 2.5.2): a
#      MASS-fraction-weighted average of Z_i/A_i, used elsewhere in this
#      project (analytic_solution.ipynb, analyze_dedx.C) for the effective
#      Z/A of water in the Bethe stopping-power formula.
#
#   2) The Bragg-KLEEMAN range rule (this file): the empirical CSDA range
#      power law R(E) = alpha * E^p, together with a material-scaling
#      relation that transfers alpha between materials via their density
#      and an ATOM-fraction-weighted average of sqrt(A) (a DIFFERENT
#      averaging from rule 1). Source for both the power law's use for
#      protons and the material-scaling relation:
#        - Water's (alpha, p): doi:10.48550/arXiv.2011.00285, CSDA
#          approximation (alpha = 0.002777 cm/MeV^p, p = 1.723) -- same
#          constants already used by generate_energy_scan.py and
#          analyze_dedx.C for the thin-slab rule.
#        - The scaling relation itself, R = R_ref * (rho_ref/rho) *
#          sqrt(A)/sqrt(A_ref), and the sqrt(A) additivity rule for
#          compounds/mixtures: W.H. Bragg, R. Kleeman, "On the alpha
#          particles of radium, and their loss of range in passing
#          through various atoms and molecules", Philosophical Magazine,
#          Series 6, 10(57), 318-340 (1905).
#
# CAVEAT (documented, not fixed): the 1905 rule was derived for alpha
# particles and is most accurate for light-to-medium-Z absorbers; applying
# it here to protons (following arXiv:2011.00285's own use of the power
# law for protons in water) and extrapolating to a high-Z material like
# lead carries a larger, undocumented systematic uncertainty than the
# water-only fit. Treat the aluminium/copper/lead constants below as a
# first-order estimate, not a substitute for a material-specific fit to
# PSTAR/CSDA data, until cross-checked.
#
# WATER'S (alpha, p) is the project's anchor point; this module does NOT
# re-derive it, only scales outward from it to the other four materials.
#
# THE SCALING RULE applies to the WHOLE range curve, not just one energy,
# because the scaling factor below is energy-independent: R = alpha*E^p
# transforms into another EXACT power law with the SAME exponent p:
#
#     R_target(E) = R_ref(E) * (rho_ref / rho_target)
#                    * sqrt(A_target) / sqrt(A_ref)
#
# with sqrt(A) for a compound/mixture given by the original Bragg-Kleeman
# additivity rule (ATOM/mole fraction weighted -- NOT mass fraction):
#
#     sqrt(A_eff) = sum_i n_i * sqrt(A_i),   n_i = atom fraction of i
#
# which reduces to sqrt(A) itself for a pure element.
#
# WHY BRIDGE THROUGH AIR: water's constants are scaled to air first
# (mirroring the rule's original 1905 form, historically expressed
# relative to air) and then from air to aluminium/copper/lead.
# Algebraically the air step cancels out (water -> X directly gives the
# same alpha_X as water -> air -> X), but routing through air lets each
# step be checked against an independently confirmed reference (air's
# density and sqrt(A), see verify_air_reference_values() below) before
# it feeds the next material.
# ============================================================================
import math

import numpy as np

# --- Standard atomic weights (IUPAC/CIAAW 2021, https://iupac.qmul.ac.uk/AtWt/),
# matching the values already used in analytic_solution.ipynb / analyze_dedx.C:
#   H  = 1.0080(2) g/mol
#   C  = 12.011    g/mol  (air trace constituent, CO2)
#   N  = 14.007    g/mol
#   O  = 15.999(1) g/mol
#   Ar = 39.948    g/mol  (air trace constituent)
#   Al = 26.981 5384(3) g/mol
#   Cu = 63.546(3) g/mol
#   Pb = 207.2(1.1) g/mol  -- large uncertainty is genuine natural isotopic
#                             variability (IUPAC 2021 interval-covering value)
ATOMIC_WEIGHT_G_MOL = {
    "H": 1.0080,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "Ar": 39.948,
    "Al": 26.9815384,
    "Cu": 63.546,
    "Pb": 207.2,
}

# --- Air composition: G4_AIR / NIST mass fractions (Z, w):
#   C(6) 0.000124, N(7) 0.755268, O(8) 0.231781, Ar(18) 0.012827
AIR_MASS_FRACTIONS = {"C": 0.000124, "N": 0.755268, "O": 0.231781, "Ar": 0.012827}

# Air density: standard dry-air value at 0 C (273.15 K), 1 atm (101325 Pa)
# -- the historical reference condition of the original Bragg-Kleeman rule
# (and of the classic alpha-particle range tables it was built for), NOT
# Geant4's default "G4_AIR" / "AIR, DRY (NEAR SEA LEVEL)" NIST material,
# which is quoted at 20 C (rho ~ 1.2048e-3 g/cm3). Cross-checked here via
# the ideal gas law with the standard molar mass of dry air
# (M = 28.9647 g/mol, US Standard Atmosphere 1976) and CODATA 2022
# R = 8.314462618 J/(mol K): rho = PM/(RT) = 1.2923e-3 g/cm3, matching the
# widely tabulated value (e.g. CRC Handbook of Chemistry and Physics) to
# better than 0.06%.
AIR_DENSITY_G_CM3 = 1.2929e-3


def _air_atom_fractions() -> dict:
    """Atom (mole) fractions of air from its NIST mass-fraction
    composition: n_i ~ w_i / A_i, normalized to sum to 1 -- the standard
    mass-fraction -> mole-fraction conversion (needed because the
    Bragg-Kleeman sqrt(A) additivity rule below is atom-count-weighted,
    while G4_AIR's composition is tabulated by mass).
    """
    n_raw = {el: w / ATOMIC_WEIGHT_G_MOL[el] for el, w in AIR_MASS_FRACTIONS.items()}
    total = sum(n_raw.values())
    return {el: n / total for el, n in n_raw.items()}


# --- Materials: density [g/cm3] (PubChem Periodic Table,
# https://pubchem.ncbi.nlm.nih.gov/ptable/density/, except air -- see
# above) and ATOM composition (element -> atom count for a stoichiometric
# compound, or atom fraction for a mixture) ----------------------------------
MATERIALS = {
    "water": {"density_g_cm3": 1.0, "atoms": {"H": 2, "O": 1}},  # H2O, exact stoichiometry
    "air": {"density_g_cm3": AIR_DENSITY_G_CM3, "atoms": _air_atom_fractions()},
    "aluminium": {"density_g_cm3": 2.70, "atoms": {"Al": 1}},
    "copper": {"density_g_cm3": 8.96, "atoms": {"Cu": 1}},
    "lead": {"density_g_cm3": 11.34, "atoms": {"Pb": 1}},
}

# Geant4 NIST material names, for use by generate_energy_scan.py.
G4_MATERIAL_NAME = {
    "water": "G4_WATER",
    "air": "G4_AIR",
    "aluminium": "G4_Al",
    "copper": "G4_Cu",
    "lead": "G4_Pb",
}

# Water's Bragg-Kleeman constants: the project's sourced anchor point.
# Source: doi:10.48550/arXiv.2011.00285, CSDA approximation.
WATER_BK_ALPHA_CM = 0.002777  # cm / MeV^p
WATER_BK_P = 1.723


def sqrt_A_effective(material_key: str) -> float:
    """
    Effective sqrt(A) of `material_key` for the Bragg-Kleeman range-scaling
    rule:

        sqrt(A_eff) = sum_i n_i * sqrt(A_i)

    with n_i the ATOM (mole) fraction of element i -- the original
    Bragg-Kleeman additivity rule for compounds/mixtures (Bragg & Kleeman,
    1905). For a pure element this reduces to sqrt(A) of that element.
    """
    atoms = MATERIALS[material_key]["atoms"]
    total = sum(atoms.values())
    return sum((n / total) * math.sqrt(ATOMIC_WEIGHT_G_MOL[el])
                for el, n in atoms.items())


def scale_bragg_kleeman_alpha(alpha_ref: float, ref_key: str, target_key: str) -> float:
    """
    Scale a Bragg-Kleeman alpha constant (R = alpha * E^p) from `ref_key`
    to `target_key`, assuming the same exponent p for both -- valid
    because the scaling factor below is energy-independent, so p
    transfers exactly:

        alpha_target = alpha_ref * (rho_ref / rho_target)
                        * sqrt(A_target) / sqrt(A_ref)
    """
    rho_ref = MATERIALS[ref_key]["density_g_cm3"]
    rho_tgt = MATERIALS[target_key]["density_g_cm3"]
    sqrtA_ref = sqrt_A_effective(ref_key)
    sqrtA_tgt = sqrt_A_effective(target_key)
    return alpha_ref * (rho_ref / rho_tgt) * (sqrtA_tgt / sqrtA_ref)


def fit_bragg_kleeman(energies_MeV, ranges_cm):
    """
    Recover (alpha, p) of R = alpha * E^p from tabulated (E, R) points via
    a log-log linear regression: ln(R) = ln(alpha) + p * ln(E). This is
    the "interpolate the constants from the computed range" step used by
    derive_all_constants() for aluminium/copper/lead below. Unlike the
    closed-form algebraic scaling used for air, this generalizes directly
    to fitting a range curve that is NOT an exact power law (e.g. actual
    simulated CSDA ranges from a future Geant4 run in each material).
    """
    log_E = np.log(np.asarray(energies_MeV, dtype=float))
    log_R = np.log(np.asarray(ranges_cm, dtype=float))
    p, log_alpha = np.polyfit(log_E, log_R, 1)
    return math.exp(log_alpha), p


def derive_all_constants(energies_MeV) -> dict:
    """
    Derive Bragg-Kleeman (alpha, p) for every material in MATERIALS from
    water's constants (arXiv:2011.00285), bridging through air:

      1) alpha_air is obtained ALGEBRAICALLY from alpha_water
         (scale_bragg_kleeman_alpha): exact, no fitting needed, since the
         scaling factor is energy-independent.
      2) R_air(E) is evaluated on `energies_MeV` with (alpha_air, p_water).
      3) For aluminium/copper/lead, R_target(E_i) is computed point by
         point from R_air(E_i) via the SAME scaling formula (forward,
         air -> target), then (alpha_target, p_target) is recovered by
         fit_bragg_kleeman() -- the "interpolate the constants from the
         computed range" step. p_target should equal p_water to
         numerical precision; this doubles as a CONSISTENCY CHECK on the
         whole pipeline; a large deviation would signal a bug, not a
         physical effect (the underlying relation is an exact power law).

    Returns {material_key: (alpha_cm, p)}. `energies_MeV` only needs to
    span a reasonably wide range with enough points for a stable fit --
    it does NOT need to match the actual simulation sweep grid.
    """
    energies_MeV = np.asarray(energies_MeV, dtype=float)
    alpha_air = scale_bragg_kleeman_alpha(WATER_BK_ALPHA_CM, "water", "air")
    R_air = alpha_air * energies_MeV ** WATER_BK_P

    constants = {
        "water": (WATER_BK_ALPHA_CM, WATER_BK_P),
        "air": (alpha_air, WATER_BK_P),
    }
    rho_air = MATERIALS["air"]["density_g_cm3"]
    sqrtA_air = sqrt_A_effective("air")
    for key in ("aluminium", "copper", "lead"):
        rho = MATERIALS[key]["density_g_cm3"]
        sqrtA = sqrt_A_effective(key)
        R_target = R_air * (rho_air / rho) * (sqrtA / sqrtA_air)
        alpha_fit, p_fit = fit_bragg_kleeman(energies_MeV, R_target)
        if abs(p_fit - WATER_BK_P) > 1e-6:
            print(f"warning: bragg_kleeman_materials: fitted p for '{key}' "
                  f"= {p_fit:.6f}, expected {WATER_BK_P:.6f} -- check the "
                  "energy grid used for the fit")
        # Keep p exactly equal to water's p (see docstring above): the fit
        # is a consistency check, not an independently meaningful value.
        constants[key] = (alpha_fit, WATER_BK_P)
    return constants


def verify_air_reference_values() -> float:
    """
    Sanity check against the values independently confirmed for this
    project from the G4_AIR NIST composition (Z, mass fraction): C(6)
    0.000124, N(7) 0.755268, O(8) 0.231781, Ar(18) 0.012827 -- density
    1.2929e-3 g/cm3, sqrt(A_air) = 3.81. NOTE: 3.81 is sqrt(A_air) itself
    (what the Bragg-Kleeman formula actually needs), not A_air -- squaring
    it gives the more familiar "mean atomic weight per atom of air"
    figure (~14.5 g/mol, consistent with air being ~78% N2/22% O2 by
    atom count, whose per-atom rather than per-molecule mean is ~half the
    ~29 g/mol mean molecular weight of air).
    """
    sqrtA_air = sqrt_A_effective("air")
    rho_air = MATERIALS["air"]["density_g_cm3"]
    assert abs(rho_air - 1.2929e-3) < 1e-9, (
        f"air density = {rho_air:.6e} g/cm3, expected 1.2929e-3 g/cm3")
    assert abs(sqrtA_air - 3.81) < 0.005, (
        f"sqrt(A_air) = {sqrtA_air:.4f}, expected ~3.81")
    return sqrtA_air


if __name__ == "__main__":
    sqrtA_air = verify_air_reference_values()
    print(f"air: rho = {MATERIALS['air']['density_g_cm3']:.6e} g/cm3 "
          f"(expected 1.2929e-03), sqrt(A_air) = {sqrtA_air:.4f} "
          "(expected ~3.81)  -- OK\n")

    # Demo: derive and print the full constants table. The fit grid here
    # is independent of any actual simulation sweep grid (see docstring).
    demo_energies = np.geomspace(3.0, 1000.0, 200)
    constants = derive_all_constants(demo_energies)
    print(f"{'material':<12}{'alpha [cm/MeV^p]':>20}{'p':>10}"
          f"{'rho [g/cm3]':>14}{'sqrt(A)':>10}")
    for key, (alpha, p) in constants.items():
        rho = MATERIALS[key]["density_g_cm3"]
        sqrtA = sqrt_A_effective(key)
        print(f"{key:<12}{alpha:>20.6e}{p:>10.4f}{rho:>14.5g}{sqrtA:>10.4f}")