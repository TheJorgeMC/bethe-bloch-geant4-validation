// ============================================================================
// analyze_dedx.C — ROOT analysis macro for the slab Bethe-Bloch validation
//
// Reads the merged CSV output of the Geant4 application (one
// "<base>_nt_slab.csv" ntuple per run and one "<base>_h1_DepthEdep.csv"
// depth-dose histogram) and produces:
//
//   1) Per-run statistics: mean/rms energy deposit (straggling), restricted
//      and total dE/dx estimators, per-event energy-balance residual, and
//      Gaussian + Landau fits of the deposit distribution.
//   2) The dE/dx vs energy curve of the full sweep, compared against the
//      PURE RELATIVISTIC BETHE formula (NO corrections) evaluated in this
//      macro — see the note below.
//   3) The Bragg / depth-dose curve from the DepthEdep histogram.
//   4) A 3D straggling surface: the per-event deposit distribution
//      (normalized to its per-run mean) for every energy of the sweep.
//   5) A cross-material error-summary table (material, actual energy range
//      analyzed, mean |rel. error| for sim-vs-Bethe, sim-vs-PSTAR and
//      Bethe-vs-PSTAR) — see ERROR SUMMARY TABLE below.
//
// UNCERTAINTY PROPAGATION (analytic side): the Bethe reference itself
// carries an uncertainty, propagated through {K, Z/A, m_e, I} (+ rho for
// the linear form), using the SAME standard uncorrelated-propagation
// formula and the SAME source values (ICRU 90 for I, IUPAC/CIAAW for the
// atomic weights behind Z/A, CODATA 2022 for K and m_e, an assumed 0.5%
// on rho) as propagate_mass_stopping_power_uncertainty() /
// propagate_linear_stopping_power_uncertainty() in analytic_solution.ipynb
// — see BetheNoCorr_*_Uncertainty_* below. dedx_summary.csv combines this
// with the simulation's own statistical error (dedxTotalErr) to report a
// propagated uncertainty on rel_err_bethe_vs_total_pct, not just the
// central value.
//
// ANALYTIC REFERENCE (consistent with analytic_solution.ipynb and the
// "Relativistic Bethe" section of Bethe_Full_Derivation.pdf):
//
//   S(T) = K z^2 (Z/A) (1/beta^2) [ ln(2 m_e c^2 beta^2 gamma^2 / I) - beta^2 ]
//
// deliberately EXCLUDING shell (C/Z), Barkas (z^3), Bloch and density
// (Sternheimer delta) corrections, and using the heavy-projectile
// approximation Tmax ~ 2 m_e c^2 beta^2 gamma^2 inside the logarithm.
// Constants are CODATA 2022 / PDG. Z/A, I and rho are MATERIAL-DEPENDENT
// (see MATERIAL AUTO-DETECTION below); water's Z/A follows the Bragg
// additivity rule with NIST STAR mass fractions, I = 78(2) eV (ICRU 90).
// EXPECTED BEHAVIOR: this curve systematically UNDERESTIMATES both NIST
// PSTAR and the Geant4 simulation, increasingly below ~10 MeV — that is a
// documented consequence of the project scope, NOT a bug.
//
// MATERIAL AUTO-DETECTION: this macro supports water, aluminium, copper
// and lead (I, Z/A, rho and their uncertainties for each — see
// kMaterials below — taken from analytic_solution.ipynb, same sourcing:
// ICRU 90/49 for I, IUPAC/CIAAW for the atomic weights behind Z/A,
// PubChem for rho). AnalyzeScan()/analyze_dedx() detect which material a
// run is by scanning the data directory for the material tag embedded in
// the file names by generate_energy_scan.py's MATERIAL_KEY (see FILE-NAME
// TAGS below); AnalyzeRun() detects it directly from the given file name.
// Air is NOT included: it has no sourced I value in this project (it is
// only a bridging material for the Bragg-Kleeman range-scaling rule, not
// a dE/dx validation target) — a directory of air data falls back to
// water with a warning.
//
// FILE-NAME TAGS: output files are named
// "dedx_<materialTag>_<Etag>MeV_<cutTag>", e.g.
// "dedx_water_0p001MeV_cut0p01mm_nt_slab.csv" — 'p'-encoded because
// G4AnalysisManager treats everything after the last '.' in a file name
// as a file-type extension. The EnergyTag() helper below must stay
// consistent with etag() in generate_energy_scan.py, and the material
// tags in kMaterials (below) with bk.G4_MATERIAL_NAME's keys in
// bragg_kleeman_materials.py.
//
// NIST PSTAR COMPARISON: in addition to the analytic Bethe reference, the
// scan is compared against the tabulated NIST PSTAR "total" stopping power
// (electronic + nuclear, stopping_power_total_MeV_cm2_g column) — the same
// quantity as the simulation's total (~unrestricted) estimator and the
// benchmark already used in analytic_solution.ipynb (100 MeV water: PSTAR
// total = 7.289 MeV cm2/g). Files are expected as
// "<pstarDir>/pstar_<materialPstarName>.csv" with columns
// material,matno,energy_MeV,stopping_power_electronic_MeV_cm2_g,
// stopping_power_nuclear_MeV_cm2_g,stopping_power_total_MeV_cm2_g,
// csda_range_g_cm2,projected_range_g_cm2,detour_factor (the NIST PSTAR
// preferred-number energy grid, same as energy_grid.csv). Note the
// materialPstarName spelling mismatch with the G4-sweep tag for aluminium:
// PSTAR uses American "aluminum" (see kMaterials' pstarStem field below).
// The PSTAR data directory is configurable — see cfg::kPstarDir below, and
// the optional pstarDir argument of AnalyzeScan()/analyze_dedx(). Missing
// or out-of-range PSTAR data degrades gracefully: the comparison is simply
// left out (dedx_summary.csv fields empty), not a fatal error. In addition
// to the two simulation comparisons (vs Bethe, vs PSTAR), dedx_summary.csv
// also carries a THEORY-VS-NIST comparison (rel_err_bethe_vs_pstar_pct/
// _err) that needs no simulation data at all — see BetheNoCorr_Mass_MeVcm2_g
// vs PSTAR directly, gated only on PSTAR coverage (not on n_exit).
//
// ENERGY RANGE: cfg::kEMinMeV / cfg::kEMaxMeV bound which grid points are
// analyzed (material detection, the sweep loop, dedx_summary.csv rows, and
// the error-summary table's reported range) WITHOUT touching kEnergiesMeV/
// energy_grid.csv itself — e.g. set kEMaxMeV = 300.0 to restrict everything
// to the clinical 3-300 MeV band. kEMaxMeV defaults to effectively "no
// upper limit" so existing behavior is unchanged unless edited.
//
// ERROR SUMMARY TABLE: every AnalyzeScan()/analyze_dedx() call appends one
// row (material, ACTUAL energy range analyzed, mean |rel. error| for each
// of the three comparisons above) to a session-wide table, printed again in
// full (and (re)written to error_summary.csv) at the end of every call —
// see g_errorSummary/PrintErrorSummaryTable() below. Running it once per
// material in the same ROOT session (see Usage) therefore builds up ONE
// consolidated table across all of them, instead of having to compare
// several dedx_summary.csv files by hand.
//
// Usage (from the directory containing the CSV files):
//   root -l -b -q analyze_dedx.C                      # full sweep + plots
//   root -l 'analyze_dedx.C("path/to/data")'          # data in another dir
//   root -l 'analyze_dedx.C("path/to/data", "path/to/nist_data")'
//                                                       # + custom PSTAR dir
//   root -l
//     .L analyze_dedx.C
//     AnalyzeRun("dedx_water_150MeV_cut0p01mm_nt_slab.csv", 150.0);  // single run
//     PlotBragg("dedx_water_150MeV_cut0p01mm_h1_DepthEdep.csv");
//     // Multiple materials in one session -> one consolidated error table:
//     AnalyzeScan("data_water"); AnalyzeScan("data_al");
//     AnalyzeScan("data_cu");    AnalyzeScan("data_pb");
//     PrintErrorSummaryTable();  // (also auto-printed after each call above)
//
// OUTPUT DIRECTORY: every file this macro writes (dedx_summary.csv, the
// PNG/PDF plots, error_summary.csv) goes into "<base>_output/", where
// <base> is the last path component of the data directory being analyzed
// (dataDir for AnalyzeScan()/analyze_dedx(); the directory of the given
// file for standalone AnalyzeRun()/PlotBragg() calls) — e.g. data in
// "path/to/data" writes outputs to "path/to/data_output/" (created
// automatically if missing, see EnsureOutputDir() below). This keeps a
// directory of raw Geant4 CSVs from being cluttered with analysis output,
// and keeps outputs from different materials/directories in separate
// folders even when run from the same working directory.
//
// The configuration block below (energy grid CSV, cut tag, thickness,
// energy range, PSTAR directory) must match the macros/data layout used to
// produce the data; the material-specific constants (Z/A, I, rho) are set
// automatically, not edited here.
// ============================================================================
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TPad.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
// ============================================================================
// Configuration — keep consistent with macros/generate_energy_scan.py and
// with analytic_solution.ipynb
// ============================================================================
static std::vector<double> readFirstColumnFromCSV(const std::string& filename)
{
  std::vector<double> result;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: could not open " << filename << "\n";
    return result;
  }
  std::string line;
  if (std::getline(file, line)) { /* skip header row */ }
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string first_column_val;
    if (std::getline(ss, first_column_val, ',')) {
      try {
        // std::stod, NOT std::stoi: the grid has decimal energies (0.001)
        // that integer parsing would silently truncate to 0.
        result.push_back(std::stod(first_column_val));
      } catch (const std::invalid_argument&) {
        // non-numeric line (e.g. stray header) — skip
      }
    }
  }
  file.close();
  return result;
}
// Energy tag for file names: %g formatting with '.' -> 'p' (0.001 -> "0p001").
// Must match etag() in generate_energy_scan.py.
static TString EnergyTag(double E_MeV)
{
  TString t = TString::Format("%g", E_MeV);
  t.ReplaceAll(".", "p");
  return t;
}
// ============================================================================
// Output directory helpers — see the OUTPUT DIRECTORY note in the header
// comment above.
// ============================================================================
// Directory portion of a file path ("a/b/c.csv" -> "a/b"; no slash -> ".").
static std::string DirOf(const std::string& path)
{
  const size_t slash = path.find_last_of('/');
  return (slash == std::string::npos) ? std::string(".") : path.substr(0, slash);
}
// Last path component of a directory, stripped of trailing slashes
// ("path/to/data/" -> "data"; "." or ".." or empty -> the fallback "data",
// since "._output"/".._output" would be a confusing folder name).
static std::string BaseNameOf(const std::string& dir)
{
  std::string d = dir;
  while (d.size() > 1 && d.back() == '/') d.pop_back();
  const size_t slash = d.find_last_of('/');
  std::string base = (slash == std::string::npos) ? d : d.substr(slash + 1);
  if (base.empty() || base == "." || base == "..") base = "data";
  return base;
}
// Computes "<basename(dataDir)>_output", creates it (recursively, no error
// if it already exists) via gSystem->mkdir(), and returns the path every
// output file should be written under.
static std::string EnsureOutputDir(const std::string& dataDir)
{
  const std::string outDir = BaseNameOf(dataDir) + "_output";
  gSystem->mkdir(outDir.c_str(), true);
  return outDir;
}
namespace cfg
{
// Minimum energy included in the analysis (MeV). Points below this are
// skipped: the uncorrected Bethe formula is outside its Born-approximation
// validity range there (see the Sommerfeld-parameter discussion in the
// derivation document) and sub-range protons stop inside the slab anyway.
const double kEMinMeV = 3.0;
// Maximum energy included in the analysis (MeV) — the upper-bound analog of
// kEMinMeV. Points above this are skipped the same way (material detection,
// the sweep loop, the summary table), letting the analysis be restricted to
// a sub-range (e.g. the clinical 3-300 MeV band) without touching
// kEnergiesMeV/energy_grid.csv itself. Effectively "no upper limit" by
// default so existing behavior is unchanged unless this is edited.
const double kEMaxMeV = 300;  //1.0e9;
// Sweep energies in MeV: NIST PSTAR preferred-number grid, same file used by
// generate_energy_scan.py (path relative to where root is executed).
const std::vector<double> kEnergiesMeV =
    readFirstColumnFromCSV("../../nist_data/energy_grid.csv");
// Cut tag in the data file names — must match CUT_TAG in
// generate_energy_scan.py ("cut0p01mm" for the 0.01 mm fine cut;
// "cut1mm" for the 1 mm default).
const char* kCutTag = "cut0p01mm";
// Ntuple file pattern: dedx_<materialTag>_<Etag>MeV_<cutTag>_nt_slab.csv
// (%s placeholders: dataDir, materialTag, EnergyTag(E), kCutTag)
const char* kNtupleFilePattern = "%s/dedx_%s_%sMeV_%s_nt_slab.csv";
// NIST PSTAR reference data directory (relative to where root is executed
// by default, same convention as the energy_grid.csv path above) — override
// via the pstarDir argument of AnalyzeScan()/analyze_dedx() if the data
// lives somewhere else. Files expected as "<kPstarDir>/pstar_<name>.csv".
const char* kPstarDir = "../../nist_data";
// SLAB THICKNESS IS NO LONGER A SINGLE NUMBER: the generator applies the
// thin-slab rule t(E) = clamp(5% x range(E), 1 um, 5 mm), with the range
// from the (material-specific) Bragg-Kleeman rule R = alpha x E^p, CSDA
// approximation, so that the fractional energy loss stays small at EVERY
// sweep energy and the total estimator (E_in - E_out)/track is defined
// across the whole validation band. The analysis needs no thickness input
// at all: every estimator divides by the per-event track length from the
// ntuple, regardless of material.
const double kThicknessMaxMM = 5.0;  // informational (upper clamp)
// --- Physical constants (CODATA 2022 / PDG, as in the notebook) -------------
const double kMp = 938.27208943;    // proton rest mass energy, MeV
const double kMe = 0.51099895069;   // electron rest mass energy, MeV
const double kK = 0.307075;         // 4 pi N_A r_e^2 m_e c^2, MeV cm2 / mol
const double kEvToMeV = 1e-6;       // eV -> MeV
const double kDeltaK = 3.0e-10;     // MeV cm2/mol, CODATA-propagated
const double kDeltaMe = 1.6e-10;    // MeV, CODATA 2022: me c^2 = ...69(16) MeV
// --- Water constituents (H2O), kept as raw element data because water is
// a COMPOUND: its Z/A and Delta(Z/A) need the two-element Bragg additivity
// rule (ICRU 37/49, Sec. 2.5.2 Eq. 2.22), unlike a pure element. NIST STAR
// mass fractions (matno 276) and IUPAC 2021 atomic weights:
//   H: w = 0.111894, A = 1.0080(2);  O: w = 0.888106, A = 15.999(1).
// Used only by the water entry of kMaterials, below.
const double kWaterWH = 0.111894, kWaterAH = 1.0080, kWaterZH = 1.0;
const double kWaterWO = 0.888106, kWaterAO = 15.999, kWaterZO = 8.0;
const double kWaterDeltaAH = 0.0002;  // IUPAC/CIAAW: A(H) = 1.0080(2) g/mol
const double kWaterDeltaAO = 0.001;   // IUPAC/CIAAW: A(O) = 15.999(1) g/mol
// --- Active material -----------------------------------------------------
// NOT hardcoded to water anymore: these are overwritten at runtime by
// SetActiveMaterial(), called by AnalyzeScan()/AnalyzeRun() once they
// detect which material a run's data belongs to (see kMaterials and the
// detection helpers below "Material properties table"). Initialized to
// water's values so anything evaluated before detection runs (there
// shouldn't be any) gets the previous water-only behavior, not garbage.
const char* kMaterialName = "WATER";
double kZoverA =
    kWaterWH * (kWaterZH / kWaterAH) + kWaterWO * (kWaterZO / kWaterAO);
double kDeltaZoverA = TMath::Sqrt(
    TMath::Sq(kWaterWH * kWaterZH / (kWaterAH * kWaterAH) * kWaterDeltaAH) +
    TMath::Sq(kWaterWO * kWaterZO / (kWaterAO * kWaterAO) * kWaterDeltaAO));
double kDensity = 1.0;       // g/cm3 (PubChem; as in the notebook)
double kDeltaDensity = kDensity * 0.005;  // 0.5% ASSUMED (no source located)
double kI_eV = 78.0;         // mean excitation energy, ICRU 90: 78(2) eV
double kDeltaI_eV = 2.0;     // ICRU 90: I(water) = 78(2) eV
// --- Straggling surface (AnalyzeScan) ---------------------------------------
// The deposit of each event is normalized to its run's mean, so runs whose
// absolute deposits differ by orders of magnitude across the sweep share
// one common, comparable y axis. 1.0 = the mean; the shape around it IS
// the straggling distribution.
// Wider y range than before: with thin slabs the high-energy runs sit deep
// in the Landau regime, whose mean-normalized delta-ray tail extends well
// above the old 1.8 limit.
const int kStragNBinsY = 80;
const double kStragYMin = 0.2;   // Edep / <Edep>
const double kStragYMax = 3.0;
}  // namespace cfg
// ============================================================================
// Material properties table: I (ICRU 90/49), Z/A (Bragg additivity rule for
// water, a compound; a pure element's Z/A is just Z/A_atomic) and their
// uncertainties, plus density (PubChem) — same values and sources as the
// MATERIALS dict in analytic_solution.ipynb. Used to auto-configure the
// cfg::k* active-material variables above for whichever material's data is
// found in a given run (SetActiveMaterial() below). "tag" is the file-name
// token generate_energy_scan.py's MATERIAL_KEY embeds in output file names
// (dedx_<tag>_<Etag>MeV_<cutTag>...).
//
// AIR IS DELIBERATELY NOT INCLUDED: this project has no sourced mean
// excitation energy I for air (analytic_solution.ipynb only tabulates I for
// water/aluminium/copper/lead), and air is not a target of the dE/dx
// validation — it is only a bridging material for the Bragg-Kleeman
// range-scaling rule in generate_energy_scan.py (see
// bragg_kleeman_materials.py). A directory of MATERIAL_KEY="air" data will
// not match any tag below and falls back to water, with a warning (see
// DetectMaterialIndexInDir()).
// ============================================================================
struct MaterialProps
{
  const char* name;              // display name (dedx_summary.csv column)
  const char* tag;                // file-name tag (matches generate_energy_scan.py)
  const char* pstarStem;          // NIST PSTAR file stem: pstar_<pstarStem>.csv
  double density, deltaDensity;   // g/cm3
  double ZoverA, deltaZoverA;     // mol/g
  double I_eV, deltaI_eV;         // eV
};
// Z/A and its uncertainty for a PURE ELEMENT: Z/A = Z/A_atomic,
// Delta(Z/A) = (Z/A)/A_atomic * Delta(A_atomic) — same formula as the
// non-compound branch of _delta_ZA() in analytic_solution.ipynb.
static MaterialProps ElementMaterial(const char* name, const char* tag,
                                     const char* pstarStem,
                                     double Z, double A, double deltaA,
                                     double density, double deltaDensity,
                                     double I_eV, double deltaI_eV)
{
  const double ZoverA = Z / A;
  const double deltaZoverA = (ZoverA / A) * deltaA;
  return MaterialProps{name, tag, pstarStem, density, deltaDensity,
                       ZoverA, deltaZoverA, I_eV, deltaI_eV};
}
// Water needs the two-element (H, O) Bragg additivity rule instead of the
// single-element formula above; reuses the raw constituent data already
// declared in cfg (kWaterWH, kWaterAH, ... kWaterDeltaAO).
static MaterialProps WaterMaterial()
{
  return MaterialProps{"WATER", "water", "water", 1.0, 1.0 * 0.005,
                       cfg::kZoverA, cfg::kDeltaZoverA, 78.0, 2.0};
}
// Sourcing for aluminium/copper/lead (identical values to
// analytic_solution.ipynb's MATERIALS dict):
//   I: ICRU Report 49 (1993), Table 2.8, Sec. 2.5.1.
//   A (standard atomic weight): IUPAC/CIAAW 2021, Table 2.
//   rho: PubChem Periodic Table (pubchem.ncbi.nlm.nih.gov/ptable/density);
//        0.5% ASSUMED uncertainty (no source located for solid densities).
// NOTE the pstarStem spelling: NIST PSTAR's own downloaded file names use
// American English ("aluminum"), while this project's own G4-sweep file-tag
// convention (the "tag" field, matching bragg_kleeman_materials.py) uses
// British English ("aluminium") — two independent, deliberately-kept
// spellings, not a typo in either place.
static const std::vector<MaterialProps> kMaterials = {
  WaterMaterial(),
  // Aluminium: Z=13, A = 26.9815384(3) g/mol, I = 166(2) eV, rho = 2.70 g/cm3.
  ElementMaterial("ALUMINIUM", "aluminium", "aluminum", 13.0, 26.9815384,
                  0.0000003, 2.70, 2.70 * 0.005, 166.0, 2.0),
  // Copper: Z=29, A = 63.546(3) g/mol, I = 322(10) eV, rho = 8.96 g/cm3.
  ElementMaterial("COPPER", "copper", "copper", 29.0, 63.546, 0.003,
                  8.96, 8.96 * 0.005, 322.0, 10.0),
  // Lead: Z=82, A = 207.2(1.1) g/mol (large uncertainty is genuine natural
  // isotopic variability, IUPAC's own interval-covering value, not a typo),
  // I = 823(30) eV, rho = 11.34 g/cm3.
  ElementMaterial("LEAD", "lead", "lead", 82.0, 207.2, 1.1,
                  11.34, 11.34 * 0.005, 823.0, 30.0),
};
// Overwrites the cfg::k* active-material variables from kMaterials[index].
void SetActiveMaterial(int index)
{
  const MaterialProps& m = kMaterials[index];
  cfg::kMaterialName = m.name;
  cfg::kDensity = m.density;
  cfg::kDeltaDensity = m.deltaDensity;
  cfg::kZoverA = m.ZoverA;
  cfg::kDeltaZoverA = m.deltaZoverA;
  cfg::kI_eV = m.I_eV;
  cfg::kDeltaI_eV = m.deltaI_eV;
}
int MaterialIndexFromTag(const std::string& tag)
{
  for (size_t m = 0; m < kMaterials.size(); ++m)
    if (tag == kMaterials[m].tag) return (int)m;
  return -1;
}
// Extracts the material tag: the token between "dedx_" and the next '_' in
// a ntuple file's base name, matching generate_energy_scan.py's convention
// dedx_<materialTag>_<Etag>MeV_<cutTag>_nt_slab.csv. Returns "" if the file
// name doesn't start with "dedx_" or has no further '_' (e.g. a pre-
// material-tag legacy file name).
std::string MaterialTagFromFilename(const std::string& path)
{
  const size_t slash = path.find_last_of('/');
  const std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
  const std::string prefix = "dedx_";
  if (base.rfind(prefix, 0) != 0) return "";
  const size_t start = prefix.size();
  const size_t us = base.find('_', start);
  if (us == std::string::npos) return "";
  return base.substr(start, us - start);
}
// Scans `dataDir` for which known material's tagged files are present (each
// candidate tag is tried against a few grid energies), used by
// AnalyzeScan()/analyze_dedx() to auto-configure the analysis without the
// caller having to specify a material by hand. Falls back to WATER (index
// 0) with a warning if nothing matches (legacy untagged file names, or
// MATERIAL_KEY = "air", which has no entry in kMaterials — see the note
// above kMaterials).
int DetectMaterialIndexInDir(const char* dataDir)
{
  std::vector<int> found;
  for (size_t m = 0; m < kMaterials.size(); ++m) {
    for (double E : cfg::kEnergiesMeV) {
      if (E < cfg::kEMinMeV || E > cfg::kEMaxMeV) continue;
      TString f = TString::Format(cfg::kNtupleFilePattern, dataDir,
                                  kMaterials[m].tag, EnergyTag(E).Data(),
                                  cfg::kCutTag);
      std::ifstream test(f.Data());
      if (test) { found.push_back((int)m); break; }
    }
  }
  if (found.empty()) {
    printf("[warning] DetectMaterialIndexInDir: no known-material data "
           "files found in '%s' (tried: water, aluminium, copper, lead) — "
           "defaulting to WATER. Check that generate_energy_scan.py's "
           "MATERIAL_KEY-tagged file names match cfg::kNtupleFilePattern.\n",
           dataDir);
    return 0;
  }
  if (found.size() > 1) {
    printf("[warning] DetectMaterialIndexInDir: data files for %zu "
           "different materials found in '%s' — using '%s'. Keep one "
           "material per directory for a clean automatic analysis.\n",
           found.size(), dataDir, kMaterials[found[0]].name);
  }
  return found[0];
}
// ============================================================================
// NIST PSTAR reference data — loader and log-log interpolation
// ============================================================================
// One tabulated (energy, total mass stopping power) point from a
// pstar_<name>.csv file (see the NIST PSTAR COMPARISON note in the header
// comment for the full column layout; only energy_MeV and
// stopping_power_total_MeV_cm2_g are needed here).
struct PstarPoint
{
  double E_MeV;
  double S_total_MeVcm2_g;
};
// Loads "<pstarDir>/pstar_<pstarStem>.csv". Proper comma-split parsing (not
// sscanf with "%lf"): the first two columns (material, matno) are strings,
// not numbers. Returns an empty vector (with a warning) if the file cannot
// be opened — callers must treat that as "PSTAR comparison unavailable for
// this material/directory", not a fatal error.
std::vector<PstarPoint> LoadPstarData(const std::string& pstarDir,
                                      const std::string& pstarStem)
{
  std::vector<PstarPoint> pts;
  const std::string path = pstarDir + "/pstar_" + pstarStem + ".csv";
  std::ifstream in(path);
  if (!in) {
    printf("[warning] LoadPstarData: cannot open '%s' — NIST PSTAR "
           "comparison will be skipped for this run. Check the pstarDir "
           "argument / cfg::kPstarDir (currently a path relative to where "
           "root is executed).\n", path.c_str());
    return pts;
  }
  std::string line;
  if (std::getline(in, line)) { /* skip header row */ }
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> cols;
    while (std::getline(ss, field, ',')) cols.push_back(field);
    // material,matno,energy_MeV,stopping_power_electronic_MeV_cm2_g,
    // stopping_power_nuclear_MeV_cm2_g,stopping_power_total_MeV_cm2_g,...
    if (cols.size() < 6) continue;
    try {
      PstarPoint p;
      p.E_MeV = std::stod(cols[2]);
      p.S_total_MeVcm2_g = std::stod(cols[5]);
      pts.push_back(p);
    } catch (const std::invalid_argument&) {
      continue;  // stray non-numeric line — skip
    }
  }
  printf("Loaded %zu NIST PSTAR points from '%s'\n", pts.size(), path.c_str());
  return pts;
}
// Log-log linear interpolation of the PSTAR total mass stopping power at
// E_MeV: ln(S) is close to linear in ln(E) locally (the underlying physics
// is close to a power law over any short interval), which is the standard
// interpolation choice for stopping-power tables spanning decades in both
// axes. Given the sweep energies already coincide with PSTAR's own grid
// (both come from the same energy_grid.csv), this is mostly a safety net
// for --log N sweeps that fall between grid points, not a big extrapolation.
// Returns false (S_out untouched) if pts has fewer than 2 entries or E_MeV
// falls outside PSTAR's own tabulated range — deliberately NOT extrapolated
// beyond NIST's own data.
bool InterpolatePstarLogLog(const std::vector<PstarPoint>& pts, double E_MeV,
                            double& S_out)
{
  if (pts.size() < 2) return false;
  if (E_MeV < pts.front().E_MeV || E_MeV > pts.back().E_MeV) return false;
  size_t hi = 0;
  while (hi < pts.size() && pts[hi].E_MeV < E_MeV) ++hi;
  if (hi == 0) { S_out = pts.front().S_total_MeVcm2_g; return true; }
  if (hi >= pts.size()) { S_out = pts.back().S_total_MeVcm2_g; return true; }
  const PstarPoint& lo = pts[hi - 1];
  const PstarPoint& up = pts[hi];
  if (E_MeV == lo.E_MeV) { S_out = lo.S_total_MeVcm2_g; return true; }
  const double lnE = std::log(E_MeV);
  const double lnElo = std::log(lo.E_MeV), lnEup = std::log(up.E_MeV);
  const double lnSlo = std::log(lo.S_total_MeVcm2_g);
  const double lnSup = std::log(up.S_total_MeVcm2_g);
  const double t = (lnE - lnElo) / (lnEup - lnElo);
  S_out = std::exp(lnSlo + t * (lnSup - lnSlo));
  return true;
}
// ============================================================================
// Analytical stopping power: pure relativistic Bethe, NO corrections
// (exactly the mass_stopping_power() of analytic_solution.ipynb)
// ============================================================================
// Mass stopping power, MeV cm2/g — the notebook's mass_stopping_power().
double BetheNoCorr_Mass_MeVcm2_g(double Ekin_MeV)
{
  // S(T) = K z^2 (Z/A) (1/beta^2) [ ln(2 m_e c^2 beta^2 gamma^2 / I) - beta^2 ]
  //
  // Deliberately WITHOUT shell, Barkas, Bloch or density-effect corrections,
  // and with the heavy-projectile approximation Tmax ~ 2 m_e c^2 b^2 g^2
  // absorbed into the single logarithm (see "Distribution of the
  // relativistic correction and assembled result" in the derivation
  // document). It therefore systematically underestimates PSTAR/Geant4,
  // increasingly below ~10 MeV — expected, documented, not a bug.
  const double gamma = 1.0 + Ekin_MeV / cfg::kMp;
  const double beta2 = 1.0 - 1.0 / (gamma * gamma);

  const double I_MeV = cfg::kI_eV * cfg::kEvToMeV;
  const double logArg = 2.0 * cfg::kMe * beta2 * gamma * gamma / I_MeV;
  const double bracket = TMath::Log(logArg) - beta2;

  // z = 1 for protons.
  return cfg::kK * cfg::kZoverA / beta2 * bracket;
}
// Linear stopping power, MeV/cm — the notebook's linear_stopping_power():
// -dE/dx = S(T) * rho.
double BetheNoCorr_Linear_MeV_cm(double Ekin_MeV)
{
  return BetheNoCorr_Mass_MeVcm2_g(Ekin_MeV) * cfg::kDensity;
}
// ============================================================================
// Analytic Bethe uncertainty propagation — same formula, variable set and
// source values as propagate_mass_stopping_power_uncertainty() /
// propagate_linear_stopping_power_uncertainty() in analytic_solution.ipynb.
// Partial derivatives below are the closed-form equivalent of that
// notebook's sympy.diff() output for
//   S = K z^2 (Z/A) (1/beta^2) [ln(2 me beta^2 gamma^2 / I) - beta^2]
// (z = 1 for protons): dS/dK and dS/d(Z/A) are just S/K and S/(Z/A) since S
// is linear in both; dS/dme and dS/dI come from differentiating the log
// term, ln(2 me beta^2 gamma^2) - ln(I), whose me- and I-dependence is
// purely logarithmic (d/dme = 1/me, d/dI = -1/I).
// ============================================================================
// Uncertainty on the MASS stopping power, MeV cm2/g.
double BetheNoCorr_Mass_Uncertainty_MeVcm2_g(double Ekin_MeV)
{
  const double gamma = 1.0 + Ekin_MeV / cfg::kMp;
  const double beta2 = 1.0 - 1.0 / (gamma * gamma);
  const double I_MeV = cfg::kI_eV * cfg::kEvToMeV;
  const double bracket =
      TMath::Log(2.0 * cfg::kMe * beta2 * gamma * gamma / I_MeV) - beta2;

  const double dS_dK = cfg::kZoverA / beta2 * bracket;                // = S / K
  const double dS_dZA = cfg::kK / beta2 * bracket;                    // = S / (Z/A)
  const double dS_dme = cfg::kK * cfg::kZoverA / beta2 * (1.0 / cfg::kMe);
  const double dS_dI = cfg::kK * cfg::kZoverA / beta2 * (-1.0 / I_MeV);

  const double term_K = dS_dK * cfg::kDeltaK;
  const double term_ZA = dS_dZA * cfg::kDeltaZoverA;
  const double term_me = dS_dme * cfg::kDeltaMe;
  const double term_I = dS_dI * (cfg::kDeltaI_eV * cfg::kEvToMeV);

  // Combined in quadrature (uncorrelated sources), not by direct sum.
  return TMath::Sqrt(term_K * term_K + term_ZA * term_ZA +
                      term_me * term_me + term_I * term_I);
}
// Uncertainty on the LINEAR stopping power, MeV/cm.
// Delta(S_linear) = sqrt( (rho*Delta(S_mass))^2 + (S_mass*Delta(rho))^2 ),
// rho uncorrelated with the mass-stopping-power sources (K, Z/A, me, I).
double BetheNoCorr_Linear_Uncertainty_MeV_cm(double Ekin_MeV)
{
  const double S_mass = BetheNoCorr_Mass_MeVcm2_g(Ekin_MeV);
  const double dS_mass = BetheNoCorr_Mass_Uncertainty_MeVcm2_g(Ekin_MeV);
  return TMath::Sqrt(TMath::Sq(cfg::kDensity * dS_mass) +
                      TMath::Sq(S_mass * cfg::kDeltaDensity));
}
// Propagated uncertainty on rel_err_bethe_vs_total_pct =
// 100*(S_sim - S_bethe)/S_bethe, treating the simulation's total
// (~unrestricted) estimator (statistical error S_simErr) and the analytic
// Bethe value (propagated error S_betheErr, from the functions above) as
// independent sources, combined via the standard uncorrelated formula:
//   d(rel)/dS_sim   =  100 / S_bethe
//   d(rel)/dS_bethe = -100 * S_sim / S_bethe^2
double RelErrPctUncertainty(double S_sim, double S_simErr, double S_bethe,
                            double S_betheErr)
{
  const double term_sim = (100. / S_bethe) * S_simErr;
  const double term_bethe = (100. * S_sim / (S_bethe * S_bethe)) * S_betheErr;
  return TMath::Sqrt(term_sim * term_sim + term_bethe * term_bethe);
}
// Propagated uncertainty on rel_err_pstar_vs_total_pct =
// 100*(S_sim - S_pstar)/S_pstar, STATISTICAL-ERROR-ONLY: unlike the analytic
// Bethe reference (RelErrPctUncertainty above, which propagates {K, Z/A,
// m_e, I}), NIST PSTAR publishes no uncertainty on its tabulated stopping
// powers, so only the simulation's own statistical error (S_simErr) is
// propagated: d(rel)/dS_sim = 100/S_pstar, S_pstar treated as exact.
double RelErrPctUncertaintyStatOnly(double S_sim, double S_simErr, double S_ref)
{
  return (100. / S_ref) * S_simErr;
}
// ============================================================================
// CSV ntuple loader (merged file produced by the application)
// Columns: E_incident, Edep_primary, Edep_secondary, E_escaped_secondary,
//          E_exit_primary, track_length_primary   (MeV / mm)
// ============================================================================
struct RunData
{
  std::vector<double> eInc, edepPrim, edepSec, eEscaped, eExit, trackLen;
  bool ok = false;
};
RunData LoadNtuple(const char* path)
{
  RunData d;
  std::ifstream in(path);
  if (!in) {
    printf("  [warning] cannot open %s — skipped\n", path);
    return d;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    double v[6];
    if (std::sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf", &v[0], &v[1],
                    &v[2], &v[3], &v[4], &v[5]) != 6)
      continue;
    d.eInc.push_back(v[0]);
    d.edepPrim.push_back(v[1]);
    d.edepSec.push_back(v[2]);
    d.eEscaped.push_back(v[3]);
    d.eExit.push_back(v[4]);
    d.trackLen.push_back(v[5]);
  }
  d.ok = !d.eInc.empty();
  return d;
}
// ============================================================================
// Landau (x) Gaussian convolution ("langaus") — the standard straggling fit
// shape in detector physics, from ROOT's official langaus.C example. The
// Landau part carries the delta-ray tail (scale parameter Width = w), the
// Gaussian part the accumulated many-collision smearing (GSigma), so the
// convolution IS the practical stand-in for the Vavilov distribution at any
// thickness regime — unlike the pure Gaussian (no tail) or pure Landau
// (no finite variance), it closes at high statistics.
// Parameters: [0] Width (Landau scale w), [1] MP (most probable value,
// mode-shift corrected internally), [2] Area (normalization), [3] GSigma.
// ============================================================================
double langaufun(double* x, double* par)
{
  const double invsq2pi = 0.3989422804014;  // 1/sqrt(2 pi)
  const double mpshift = -0.22278298;       // Landau location -> true mode
  const int np = 100;                        // convolution steps
  const double sc = 5.0;                     // integration range in GSigma

  // Shift so that par[1] is the true most probable value of the Landau.
  const double mpc = par[1] - mpshift * par[0];
  const double xlow = x[0] - sc * par[3];
  const double xupp = x[0] + sc * par[3];
  const double step = (xupp - xlow) / np;

  double sum = 0.;
  for (int i = 1; i <= np / 2; ++i) {
    double xx = xlow + (i - 0.5) * step;
    sum += TMath::Landau(xx, mpc, par[0]) / par[0] *
           TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (i - 0.5) * step;
    sum += TMath::Landau(xx, mpc, par[0]) / par[0] *
           TMath::Gaus(x[0], xx, par[3]);
  }
  return par[2] * step * sum * invsq2pi / par[3];
}

// Mean and rms of a vector.
void MeanRms(const std::vector<double>& v, double& mean, double& rms)
{
  mean = rms = 0.0;
  if (v.empty()) return;
  double s = 0., s2 = 0.;
  for (double x : v) { s += x; s2 += x * x; }
  const double n = (double)v.size();
  mean = s / n;
  const double var = s2 / n - mean * mean;
  rms = (var > 0.) ? TMath::Sqrt(var) : 0.;
}
// ============================================================================
// Per-run analysis: statistics, dE/dx estimators, energy balance
// ============================================================================
// ----------------------------------------------------------------------------
// ESTIMATOR DEFINITIONS (linear stopping power, MeV/cm; ntuple lengths are
// in mm, converted in AnalyzeRunData; mass form = linear / density):
//
//  RESTRICTED  = <Edep_primary / track_length>. The primary's continuous
//    energy loss per path: excludes the energy carried away by delta rays
//    generated above the production cut. This is the Geant4 restricted
//    dE/dx with T_cut set by /physics/absorberCut — it DEPENDS on the cut
//    and is NOT the quantity the Bethe formula predicts.
//
//  TOTAL (~unrestricted) = <(E_inc - E_exit) / track_length>. ALL the
//    energy the primary lost per unit path, wherever it ended up
//    (local deposit + deltas + escaping radiation). This is the quantity
//    comparable with the unrestricted Bethe formula and with NIST PSTAR.
//    It requires the primary to EXIT the slab, which the thin-slab
//    thickness rule t(E) ~ 5% x range(E) guarantees at every sweep
//    energy; n_exit == 0 in a run now signals a broken thickness choice
//    (or E below ~0.3 MeV, where even the 1 um floor exceeds the range)
//    and is flagged, not silently zeroed.
// ----------------------------------------------------------------------------
struct RunResult
{
  double E = 0.;                    // nominal beam energy (MeV)
  double dedxRestricted = 0.;       // restricted estimator (MeV/cm)
  double dedxRestrictedErr = 0.;
  double dedxTotal = 0.;            // total (~unrestricted) estimator (MeV/cm)
  double dedxTotalErr = 0.;
  long nExit = 0;                   // events where the primary exited the slab
  double stragglingRms = 0.;        // rms of the total energy deposit
  double balanceResidual = 0.;      // <E_inc - (all channels)>  (MeV)
  long n = 0;
  bool ok = false;
};
RunResult AnalyzeRunData(const RunData& d, double Enominal)
{
  RunResult r;
  r.E = Enominal;
  if (!d.ok) return r;
  r.n = (long)d.eInc.size();
  // Per-event ratio estimators: dividing by the actual track length per
  // event (instead of the geometric thickness) removes the path-lengthening
  // bias from multiple scattering.
  std::vector<double> restricted, total, edepTotal, residual;
  restricted.reserve(r.n);
  total.reserve(r.n);
  edepTotal.reserve(r.n);
  residual.reserve(r.n);
  // Ntuple energies are in MeV and lengths in mm: the factor 10 converts
  // MeV/mm -> MeV/cm, the linear stopping power unit reported everywhere.
  const double kMmToCm = 10.0;
  for (size_t i = 0; i < d.eInc.size(); ++i) {
    if (d.trackLen[i] > 0.)
      restricted.push_back(kMmToCm * d.edepPrim[i] / d.trackLen[i]);
    if (d.eExit[i] >= 0. && d.trackLen[i] > 0.)
      total.push_back(kMmToCm * (d.eInc[i] - d.eExit[i]) / d.trackLen[i]);
    edepTotal.push_back(d.edepPrim[i] + d.edepSec[i]);
    const double exitE = (d.eExit[i] >= 0.) ? d.eExit[i] : 0.;
    residual.push_back(d.eInc[i] - (d.edepPrim[i] + d.edepSec[i] +
                                    d.eEscaped[i] + exitE));
  }
  r.nExit = (long)total.size();
  double mean, rms;
  MeanRms(restricted, mean, rms);
  r.dedxRestricted = mean;
  r.dedxRestrictedErr =
      restricted.empty() ? 0. : rms / TMath::Sqrt((double)restricted.size());
  MeanRms(total, mean, rms);
  r.dedxTotal = mean;
  r.dedxTotalErr = total.empty() ? 0. : rms / TMath::Sqrt((double)total.size());
  MeanRms(edepTotal, mean, rms);
  r.stragglingRms = rms;
  MeanRms(residual, mean, rms);
  r.balanceResidual = mean;
  r.ok = true;
  return r;
}
// Standalone single-run analysis with a canvas of the deposit distribution.
void AnalyzeRun(const char* csvFile, double Enominal)
{
  // Auto-detect the material from the file name (dedx_<materialTag>_...)
  // and configure cfg::k* accordingly, so a standalone AnalyzeRun() call
  // (e.g. from an interactive ROOT session, without going through
  // AnalyzeScan()/analyze_dedx()) still compares against the right Bethe
  // curve instead of silently defaulting to water.
  const std::string materialTag = MaterialTagFromFilename(csvFile);
  const int matIdx = materialTag.empty() ? -1 : MaterialIndexFromTag(materialTag);
  if (matIdx >= 0) {
    SetActiveMaterial(matIdx);
  } else {
    printf("[warning] AnalyzeRun: could not detect the material from '%s' "
           "(expected a 'dedx_<materialTag>_...' file name) — keeping the "
           "currently active material (%s)\n", csvFile, cfg::kMaterialName);
  }
  const std::string outDir = EnsureOutputDir(DirOf(csvFile));
  RunData d = LoadNtuple(csvFile);
  if (!d.ok) return;
  RunResult r = AnalyzeRunData(d, Enominal);
  printf("\n--- %s (material: %s) ---\n", csvFile, cfg::kMaterialName);
  printf("  events                  : %ld\n", r.n);
  printf("  exiting primaries       : %ld / %ld\n", r.nExit, r.n);
  printf("  dE/dx restricted (sim)  : %.4f +- %.4f MeV/cm  "
         "(S = %.4f MeV cm2/g)\n",
         r.dedxRestricted, r.dedxRestrictedErr,
         r.dedxRestricted / cfg::kDensity);
  printf("  dE/dx total      (sim)  : %.4f +- %.4f MeV/cm  "
         "(S = %.4f MeV cm2/g)\n",
         r.dedxTotal, r.dedxTotalErr, r.dedxTotal / cfg::kDensity);
  printf("  Bethe, no corrections   : %.4f +- %.4f MeV/cm  "
         "(S = %.4f +- %.4f MeV cm2/g)\n",
         BetheNoCorr_Linear_MeV_cm(Enominal),
         BetheNoCorr_Linear_Uncertainty_MeV_cm(Enominal),
         BetheNoCorr_Mass_MeVcm2_g(Enominal),
         BetheNoCorr_Mass_Uncertainty_MeVcm2_g(Enominal));
  printf("  straggling (rms Edep)   : %.4f MeV\n", r.stragglingRms);
  printf("  <balance residual>      : %.4g MeV\n", r.balanceResidual);
  // Distribution of the total energy deposit (straggling shape).
  double mean, rms;
  std::vector<double> edepTotal;
  edepTotal.reserve(d.eInc.size());
  for (size_t i = 0; i < d.eInc.size(); ++i)
    edepTotal.push_back(d.edepPrim[i] + d.edepSec[i]);
  MeanRms(edepTotal, mean, rms);
  const TString tag = EnergyTag(Enominal);
  const double lo = mean - 6. * rms, hi = mean + 6. * rms;
  TH1D* h = new TH1D(Form("hEdep_%s", tag.Data()),
                     Form("Energy deposit, E = %g MeV;E_{dep} in slab (MeV);events", Enominal),
                     120, lo, hi);
  for (double x : edepTotal) h->Fill(x);

  // ---------------------------------------------------------------------
  // Straggling-shape fit: Landau (x) Gaussian convolution ("langaus").
  // The energy-loss distribution in a slab interpolates between the two
  // classical limits of Landau-Vavilov theory, governed by the kappa
  // parameter (mean energy loss / Tmax): thin absorber (kappa << 1) ->
  // asymmetric Landau with a delta-ray tail; thick absorber (kappa >> 1)
  // -> Gaussian (Bohr straggling). The langaus convolution covers the
  // whole range: its Width (w) is the Landau SCALE parameter (sets the
  // peak width, FWHM ~ 4.02 w for the pure-Landau part, and the weight of
  // the high-side delta-ray tail; NOT a standard deviation — a pure Landau
  // has no finite variance), and GSigma is the Gaussian smearing from the
  // accumulated soft collisions. Pure-Gaussian and pure-Landau fits were
  // dropped: at 100k events both were rejected (chi2/ndf ~ 121 and ~ 524)
  // because neither shape is the true model and the statistical bin errors
  // shrink as 1/sqrt(N), letting the systematic mismatch dominate.
  //
  // GOODNESS OF FIT: chi2 = sum over bins of (data_i - fit_i)^2 / sigma_i^2,
  // ndf = (fitted bins) - (free parameters). A model consistent with the
  // data at the level of its statistical fluctuations gives chi2/ndf ~ 1;
  // chi2/ndf >> 1 means the SHAPE is wrong (not just noisy), chi2/ndf << 1
  // usually means overestimated bin errors. The quantitative criterion is
  // the p-value, TMath::Prob(chi2, ndf): p in ~O(0.1-1) = compatible,
  // p < ~1e-3 = shape genuinely rejected by the data.
  // ---------------------------------------------------------------------
  TF1* fLangaus = new TF1(Form("fLangaus_%s", tag.Data()), langaufun, lo, hi, 4);
  fLangaus->SetParNames("Width", "MP", "Area", "GSigma");
  // Seeds: narrow Landau + Gaussian smearing sharing the observed rms;
  // Area = histogram integral x bin width (langaufun is density-normalized).
  fLangaus->SetParameters(0.15 * rms, mean - 0.2 * rms,
                          h->Integral() * h->GetBinWidth(1), 0.7 * rms);
  fLangaus->SetParLimits(0, 1e-6, 5. * rms);   // Width > 0
  fLangaus->SetParLimits(3, 1e-6, 5. * rms);   // GSigma > 0
  h->Fit(fLangaus, "Q0R");  // Q: quiet, 0: don't auto-draw, R: fit range

  const double lgChi2 =
      (fLangaus->GetNDF() > 0)
          ? fLangaus->GetChisquare() / fLangaus->GetNDF()
          : 0.;
  const double lgProb =
      TMath::Prob(fLangaus->GetChisquare(), fLangaus->GetNDF());
  printf("  Langaus fit             : MP = %.4f +- %.4f MeV, "
         "w = %.4f +- %.4f MeV, GSigma = %.4f +- %.4f MeV, "
         "chi2/ndf = %.2f (p = %.3g)\n",
         fLangaus->GetParameter(1), fLangaus->GetParError(1),
         fLangaus->GetParameter(0), fLangaus->GetParError(0),
         fLangaus->GetParameter(3), fLangaus->GetParError(3), lgChi2, lgProb);

  TCanvas* c = new TCanvas(Form("cEdep_%s", tag.Data()), "Energy deposit", 900, 650);
  h->SetLineColor(kAzure + 2);
  h->SetLineWidth(2);
  h->SetStats(false);
  // Headroom above the peak: the top band of the pad is reserved for the
  // legend (left) and the parameter box (right), so nothing overlaps the
  // data, the curves or the markers.
  h->SetMaximum(1.65 * h->GetMaximum());
  h->Draw("hist");
  fLangaus->SetLineColor(kOrange + 7);
  fLangaus->SetLineWidth(2);
  fLangaus->SetNpx(400);  // smooth curve for the convolution
  fLangaus->Draw("same");

  // --- Central-value markers with their errors as horizontal bars --------
  // Sample mean (the estimator actually used for dE/dx) and the langaus
  // peak (maximum of the fitted convolution, with the MP fit error), each
  // drawn at its curve/histogram height with its uncertainty as an x-error
  // bar.
  const double meanErr = rms / TMath::Sqrt((double)edepTotal.size());
  TGraphErrors* gMeanPt = new TGraphErrors(1);
  gMeanPt->SetPoint(0, mean, h->GetBinContent(h->FindBin(mean)));
  gMeanPt->SetPointError(0, meanErr, 0.);
  gMeanPt->SetMarkerStyle(20);
  gMeanPt->SetMarkerSize(1.4);
  gMeanPt->SetMarkerColor(kAzure + 2);
  gMeanPt->SetLineColor(kAzure + 2);
  gMeanPt->SetLineWidth(2);
  gMeanPt->Draw("P same");

  const double peakX = fLangaus->GetMaximumX(lo, hi);
  TGraphErrors* gPeakPt = new TGraphErrors(1);
  gPeakPt->SetPoint(0, peakX, fLangaus->Eval(peakX));
  gPeakPt->SetPointError(0, fLangaus->GetParError(1), 0.);
  gPeakPt->SetMarkerStyle(22);
  gPeakPt->SetMarkerSize(1.5);
  gPeakPt->SetMarkerColor(kOrange + 7);
  gPeakPt->SetLineColor(kOrange + 7);
  gPeakPt->SetLineWidth(2);
  gPeakPt->Draw("P same");

  // --- Legend (top-left): curves + langaus markers ------------------------
  // "p" entries render the marker symbol of the object; the previous "pe"
  // option suppressed the symbols in some ROOT versions, which is why the
  // square/circle/triangle were missing from the label box.
  TLegend* legFit = new TLegend(0.12, 0.68, 0.52, 0.89);
  legFit->SetBorderSize(0);
  legFit->SetFillStyle(0);
  legFit->SetTextSize(0.030);
  legFit->AddEntry(h, "Simulation", "l");
  legFit->AddEntry(fLangaus,
                   Form("Langaus fit  (#chi^{2}/ndf = %.1f)", lgChi2), "l");
  legFit->AddEntry(gMeanPt, "sample mean", "p");
  legFit->AddEntry(gPeakPt, "langaus peak (MPV)", "p");
  legFit->Draw();

  // --- Parameter box (top-right): one line per parameter -----------------
  TPaveText* pav = new TPaveText(0.55, 0.62, 0.89, 0.89, "NDC");
  pav->SetBorderSize(0);
  pav->SetFillStyle(0);
  pav->SetTextAlign(12);  // left-adjusted
  pav->SetTextSize(0.028);
  // 4 decimals on the uncertainties: with 100k events the fit errors drop
  // below 0.001 MeV and 3 decimals would print a misleading "0.000".
  pav->AddText(Form("mean = %.4f #pm %.4f MeV", mean, meanErr));
  pav->AddText(Form("MP_{langaus} = %.4f #pm %.4f MeV",
                    fLangaus->GetParameter(1), fLangaus->GetParError(1)));
  pav->AddText(Form("w_{Landau} = %.4f #pm %.4f MeV",
                    fLangaus->GetParameter(0), fLangaus->GetParError(0)));
  pav->AddText(Form("#sigma_{Gaus} = %.4f #pm %.4f MeV",
                    fLangaus->GetParameter(3), fLangaus->GetParError(3)));
  pav->AddText(Form("p-value = %.3g", lgProb));
  pav->Draw();

  // 'p'-encoded tag in the output name too, for consistent sorting/globbing.
  c->SaveAs(Form("%s/edep_distribution_%sMeV.png", outDir.c_str(), tag.Data()));
}
// ============================================================================
// Depth-dose (Bragg) curve from the Geant4 CSV histogram file
// (g4tools wcsv format: '#' header lines with "#axis fixed N min max",
//  a column-title line, then N+2 rows "entries,Sw,Sw2,Sxw,Sx2w" where the
//  first and last rows are the under/overflow bins)
// ============================================================================
void PlotBragg(const char* h1File)
{
  const std::string outDir = EnsureOutputDir(DirOf(h1File));
  std::ifstream in(h1File);
  if (!in) {
    printf("[warning] cannot open %s — Bragg plot skipped\n", h1File);
    return;
  }
  int nbins = 0;
  double zmin = 0., zmax = 0.;
  std::vector<double> sw;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      // Parse "#axis fixed <nbins> <min> <max>"
      if (line.rfind("#axis fixed", 0) == 0) {
        std::istringstream ss(line.substr(11));
        ss >> nbins >> zmin >> zmax;
      }
      continue;
    }
    // Data rows are pure comma-separated numbers; the column-title line
    // ("entries,Sw,...") is filtered out by requiring a numeric start.
    if (!(std::isdigit((unsigned char)line[0]) || line[0] == '-' || line[0] == '.'))
      continue;
    double entries, Sw, Sw2, Sxw, Sx2w;
    if (std::sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf", &entries, &Sw, &Sw2,
                    &Sxw, &Sx2w) >= 2)
      sw.push_back(Sw);
  }
  if (nbins <= 0 || sw.empty()) {
    printf("[warning] %s: unrecognized histogram format\n", h1File);
    return;
  }
  // Strip under/overflow when present.
  int offset = 0;
  if ((int)sw.size() == nbins + 2) offset = 1;
  else if ((int)sw.size() != nbins) {
    printf("[warning] %s: %zu rows for %d bins — using min\n", h1File,
           sw.size(), nbins);
  }
  TH1D* h = new TH1D("hBragg", "Depth dose;z (mm);energy deposited (MeV)",
                     nbins, zmin, zmax);
  for (int i = 0; i < nbins && i + offset < (int)sw.size(); ++i)
    h->SetBinContent(i + 1, sw[i + offset]);
  TCanvas* c = new TCanvas("cBragg", "Depth dose", 800, 600);
  h->SetLineColor(kOrange + 7);
  h->SetLineWidth(2);
  h->SetStats(false);
  h->Draw("hist");
  const std::string braggPath = outDir + "/bragg_curve.png";
  c->SaveAs(braggPath.c_str());
  printf("Bragg curve written to %s (%d bins, z in [%g, %g] mm)\n",
         braggPath.c_str(), nbins, zmin, zmax);
}
// ============================================================================
// Cross-material error-summary table: one row per AnalyzeScan()/
// analyze_dedx() call in the current ROOT session, so running it once per
// material (water, then aluminium, then copper, then lead, ...) and calling
// PrintErrorSummaryTable() at the end (done automatically after each
// AnalyzeScan() — see below) gives a single consolidated table/CSV instead
// of having to cross-reference several dedx_summary.csv files by hand.
// ============================================================================
struct ErrorSummaryRow
{
  std::string material;
  double eMinMeV = 0., eMaxMeV = 0.;   // ACTUAL range analyzed (after
                                        // kEMinMeV/kEMaxMeV and available data)
  long nPoints = 0;                    // sweep points with readable data
  // Mean of |rel_err_*_pct| over the points where that comparison is
  // defined (see the gating conditions of each column in AnalyzeScan's CSV
  // block) — each with its own valid-point count, since the three
  // comparisons are independently gated. Naming spells out WHICH TWO
  // quantities are being compared, to avoid the ambiguity of the earlier
  // "vs total"/"vs Bethe" names (which "total" — simulation's or NIST's?):
  //   Geant4  = this project's Geant4 simulation (total/~unrestricted
  //             estimator, dedx_total_MeV_cm in dedx_summary.csv)
  //   Teoria  = the pure relativistic Bethe formula, NO corrections,
  //             evaluated analytically in this macro (bethe_no_corr_*)
  //   NIST    = the tabulated NIST PSTAR "total" stopping power
  //             (pstar_*_MeV_cm2_g in dedx_summary.csv)
  double meanAbsGeant4VsTeoriaPct = 0.;  long nGeant4VsTeoria = 0;
  double meanAbsGeant4VsNistPct = 0.;    long nGeant4VsNist = 0;
  double meanAbsTeoriaVsNistPct = 0.;    long nTeoriaVsNist = 0;
};
static std::vector<ErrorSummaryRow> g_errorSummary;
// Prints (and writes "<outDir>/error_summary.csv" with) every
// ErrorSummaryRow accumulated so far in this ROOT session. Safe to call
// more than once — it always reflects the full g_errorSummary vector, not
// just the latest run; outDir only controls WHERE the (full, cumulative)
// CSV is (re)written each time — see the OUTPUT DIRECTORY note in the
// header comment.
void PrintErrorSummaryTable(const std::string& outDir = ".")
{
  if (g_errorSummary.empty()) {
    printf("\n[error summary] no runs recorded yet — call "
           "AnalyzeScan()/analyze_dedx() first\n");
    return;
  }
  printf("\n=== Error summary across %zu analyzed run(s) ===\n"
         "    Geant4 = this project's Geant4 simulation (total estimator)\n"
         "    Teoria = pure relativistic Bethe formula (no corrections, "
         "analytic)\n"
         "    NIST   = NIST PSTAR tabulated data (total stopping power)\n\n",
         g_errorSummary.size());
  printf("%-10s %16s %7s  %-22s %-22s %-22s\n", "material", "E range (MeV)",
         "n_pts", "<|Geant4-Teoria|> %", "<|Geant4-NIST|> %",
         "<|Teoria-NIST|> %");
  const std::string csvPath = outDir + "/error_summary.csv";
  std::ofstream fout(csvPath);
  fout << "material,E_min_MeV,E_max_MeV,n_points,"
          "mean_abs_rel_err_geant4_vs_teoria_pct,n_geant4_vs_teoria,"
          "mean_abs_rel_err_geant4_vs_nist_pct,n_geant4_vs_nist,"
          "mean_abs_rel_err_teoria_vs_nist_pct,n_teoria_vs_nist\n";
  for (const ErrorSummaryRow& row : g_errorSummary) {
    char eRange[32];
    snprintf(eRange, sizeof(eRange), "%g-%g", row.eMinMeV, row.eMaxMeV);
    char sG4Teoria[40], sG4Nist[40], sTeoriaNist[40];
    if (row.nGeant4VsTeoria > 0)
      snprintf(sG4Teoria, sizeof(sG4Teoria), "%.3f (n=%ld)",
               row.meanAbsGeant4VsTeoriaPct, row.nGeant4VsTeoria);
    else
      snprintf(sG4Teoria, sizeof(sG4Teoria), "n/a");
    if (row.nGeant4VsNist > 0)
      snprintf(sG4Nist, sizeof(sG4Nist), "%.3f (n=%ld)",
               row.meanAbsGeant4VsNistPct, row.nGeant4VsNist);
    else
      snprintf(sG4Nist, sizeof(sG4Nist), "n/a");
    if (row.nTeoriaVsNist > 0)
      snprintf(sTeoriaNist, sizeof(sTeoriaNist), "%.3f (n=%ld)",
               row.meanAbsTeoriaVsNistPct, row.nTeoriaVsNist);
    else
      snprintf(sTeoriaNist, sizeof(sTeoriaNist), "n/a");
    printf("%-10s %16s %7ld  %-22s %-22s %-22s\n", row.material.c_str(),
           eRange, row.nPoints, sG4Teoria, sG4Nist, sTeoriaNist);
    fout << row.material << ',' << row.eMinMeV << ',' << row.eMaxMeV << ','
         << row.nPoints << ',';
    if (row.nGeant4VsTeoria > 0) fout << row.meanAbsGeant4VsTeoriaPct;
    fout << ',' << row.nGeant4VsTeoria << ',';
    if (row.nGeant4VsNist > 0) fout << row.meanAbsGeant4VsNistPct;
    fout << ',' << row.nGeant4VsNist << ',';
    if (row.nTeoriaVsNist > 0) fout << row.meanAbsTeoriaVsNistPct;
    fout << ',' << row.nTeoriaVsNist << '\n';
  }
  fout.close();
  printf("Written: %s (%zu material run(s) so far in this session)\n",
         csvPath.c_str(), g_errorSummary.size());
}
// ============================================================================
// Full sweep: dE/dx vs E against the uncorrected relativistic Bethe curve,
// plus the 3D straggling surface across the sweep
// ============================================================================
void AnalyzeScan(const char* dataDir = ".", const char* pstarDir = cfg::kPstarDir)
{
  gStyle->SetOptStat(0);
  if (cfg::kEnergiesMeV.empty()) {
    printf("Energy grid is empty — check the path of energy_grid.csv "
           "(read relative to the current directory)\n");
    return;
  }
  // --- Output directory: "<basename(dataDir)>_output/" (see the OUTPUT
  // DIRECTORY note in the header comment) — created up front so every
  // output below (plots, dedx_summary.csv, error_summary.csv) lands there.
  const std::string outDir = EnsureOutputDir(dataDir);
  printf("Output directory: %s\n", outDir.c_str());
  // --- Auto-detect the material and configure cfg::k* accordingly --------
  const int matIdx = DetectMaterialIndexInDir(dataDir);
  SetActiveMaterial(matIdx);
  printf("Detected material: %s (tag '%s'), I = %g(%g) eV, "
         "Z/A = %.6f(%.6f), rho = %g(%g) g/cm3\n",
         cfg::kMaterialName, kMaterials[matIdx].tag, cfg::kI_eV,
         cfg::kDeltaI_eV, cfg::kZoverA, cfg::kDeltaZoverA, cfg::kDensity,
         cfg::kDeltaDensity);
  printf("analytic Bethe uncertainty at 100 MeV (sanity check): "
         "S = %.4f +- %.4f MeV cm2/g (%.2f%% relative)\n\n",
         BetheNoCorr_Mass_MeVcm2_g(100.0),
         BetheNoCorr_Mass_Uncertainty_MeVcm2_g(100.0),
         100. * BetheNoCorr_Mass_Uncertainty_MeVcm2_g(100.0) /
             BetheNoCorr_Mass_MeVcm2_g(100.0));

  // --- NIST PSTAR reference data for the detected material ----------------
  // Degrades gracefully: if the file is missing/unreadable, havePstar stays
  // false and the whole PSTAR comparison (plot curve, ratio series,
  // dedx_summary.csv columns) is simply left out — not a fatal error.
  const std::vector<PstarPoint> pstarData =
      LoadPstarData(pstarDir, kMaterials[matIdx].pstarStem);
  const bool havePstar = pstarData.size() >= 2;
  if (havePstar) {
    printf("PSTAR comparison enabled (%zu points, E in [%g, %g] MeV)\n\n",
           pstarData.size(), pstarData.front().E_MeV, pstarData.back().E_MeV);
  } else {
    printf("PSTAR comparison DISABLED for this run (no usable data for "
           "'%s' in '%s')\n\n", kMaterials[matIdx].pstarStem, pstarDir);
  }

  // --- Pre-scan: which sweep points actually have data files? -------------
  // Needed to give the straggling TH2 one x bin per available energy.
  std::vector<double> present;
  for (double E : cfg::kEnergiesMeV) {
    if (E < cfg::kEMinMeV || E > cfg::kEMaxMeV) continue;
    TString f = TString::Format(cfg::kNtupleFilePattern, dataDir,
                                kMaterials[matIdx].tag, EnergyTag(E).Data(),
                                cfg::kCutTag);
    std::ifstream test(f.Data());
    if (test) present.push_back(E);
  }
  if (present.empty()) {
    printf("No data files found in '%s' for material '%s' — check "
           "cfg::kNtupleFilePattern and cfg::kCutTag ('%s')\n",
           dataDir, cfg::kMaterialName, cfg::kCutTag);
    return;
  }

  // --- Straggling surface: one x bin per sweep energy, y = Edep/<Edep> ----
  // Normalizing each event's deposit to its run's mean makes runs whose
  // absolute deposits differ by orders of magnitude share one comparable
  // axis: the spread/asymmetry around 1.0 IS the straggling shape, and its
  // evolution along x shows the Vavilov thick->thin transition directly.
  TH2D* hStrag = new TH2D(
      "hStrag",
      "Straggling across the energy sweep;beam energy (MeV);"
      "E_{dep}/#LTE_{dep}#GT;events",
      (int)present.size(), 0., (double)present.size(),
      cfg::kStragNBinsY, cfg::kStragYMin, cfg::kStragYMax);
  // Label only 4 representative energies (log-spaced across the sweep):
  // 100+ per-bin labels smear into an unreadable band. Every bin keeps its
  // data — unlabeled bins simply get an empty label string.
  {
    for (size_t i = 0; i < present.size(); ++i)
      hStrag->GetXaxis()->SetBinLabel((int)i + 1, "");
    const int nLab = 4;
    const int nPts = (int)present.size();
    for (int k = 0; k < nLab; ++k) {
      int idx = (nLab > 1)
                    ? (int)std::lround(k * (nPts - 1) / (double)(nLab - 1))
                    : 0;
      hStrag->GetXaxis()->SetBinLabel(idx + 1, Form("%g", present[idx]));
    }
    hStrag->GetXaxis()->SetLabelSize(0.045);
  }

  std::vector<RunResult> results;
  int ix = 0;
  for (double E : present) {
    TString f = TString::Format(cfg::kNtupleFilePattern, dataDir,
                                kMaterials[matIdx].tag, EnergyTag(E).Data(),
                                cfg::kCutTag);
    RunData d = LoadNtuple(f.Data());
    RunResult r = AnalyzeRunData(d, E);
    ++ix;
    if (!r.ok) continue;
    results.push_back(r);
    if (r.nExit > 0) {
      printf("E = %8g MeV : n = %6ld  dE/dx(restr) = %10.4f  "
             "dE/dx(total) = %10.4f  Bethe(no corr) = %10.4f +- %.4f MeV/cm\n",
             E, r.n, r.dedxRestricted, r.dedxTotal,
             BetheNoCorr_Linear_MeV_cm(E), BetheNoCorr_Linear_Uncertainty_MeV_cm(E));
    } else {
      // No exiting primaries: with the thin-slab rule this should NOT
      // happen in the validation band — flag the run for re-simulation.
      printf("E = %8g MeV : n = %6ld  dE/dx(restr) = %10.4f  "
             "dE/dx(total) =        n/a  Bethe(no corr) = %10.4f MeV/cm"
             "  [WARNING: no exiting primaries — thickness incompatible "
             "with this energy, re-simulate this point]\n",
             E, r.n, r.dedxRestricted, BetheNoCorr_Linear_MeV_cm(E));
    }
    // Fill the straggling slice for this energy.
    double meanEdep = 0., rmsEdep = 0.;
    std::vector<double> edepTotal;
    edepTotal.reserve(d.eInc.size());
    for (size_t i = 0; i < d.eInc.size(); ++i)
      edepTotal.push_back(d.edepPrim[i] + d.edepSec[i]);
    MeanRms(edepTotal, meanEdep, rmsEdep);
    if (meanEdep > 0.) {
      for (double e : edepTotal)
        hStrag->Fill(ix - 0.5, e / meanEdep);
    }
  }
  if (results.empty()) {
    printf("No readable data — nothing to analyze\n");
    return;
  }
  // --- PSTAR value at each swept energy (log-log interpolated) ------------
  std::vector<double> pstarS_MeVcm2_g(results.size(), 0.0);
  std::vector<bool> pstarOk(results.size(), false);
  if (havePstar) {
    for (size_t i = 0; i < results.size(); ++i) {
      double S = 0.;
      if (InterpolatePstarLogLog(pstarData, results[i].E, S)) {
        pstarS_MeVcm2_g[i] = S;
        pstarOk[i] = true;
      }
    }
  }
  // --- Graphs ---------------------------------------------------------------
  const int n = (int)results.size();
  TGraphErrors* gRestr = new TGraphErrors(n);
  TGraphErrors* gTotal = new TGraphErrors(n);
  TGraphErrors* gRatio = new TGraphErrors(n);
  int nRatio = 0;
  for (int i = 0; i < n; ++i) {
    const RunResult& r = results[i];
    gRestr->SetPoint(i, r.E, r.dedxRestricted);
    gRestr->SetPointError(i, 0., r.dedxRestrictedErr);
    gTotal->SetPoint(i, r.E, r.dedxTotal);
    gTotal->SetPointError(i, 0., r.dedxTotalErr);
    // With the thin-slab thickness rule the total estimator exists at every
    // sweep energy; nExit == 0 here is a data-quality flag (wrong thickness
    // for that run, or E below the ~0.3 MeV floor), not an expected gap.
    if (r.nExit > 0) {
      const double bb = BetheNoCorr_Linear_MeV_cm(r.E);
      gRatio->SetPoint(nRatio, r.E, r.dedxTotal / bb);
      gRatio->SetPointError(nRatio, 0., r.dedxTotalErr / bb);
      ++nRatio;
    }
  }
  gRatio->Set(nRatio);
  // Smooth analytical curve over the measured range.
  const double Elo = results.front().E * 0.8;
  const double Ehi = results.back().E * 1.2;
  const int nCurve = 400;
  TGraph* gBB = new TGraph(nCurve);
  for (int i = 0; i < nCurve; ++i) {
    const double E = Elo * TMath::Power(Ehi / Elo, (double)i / (nCurve - 1));
    gBB->SetPoint(i, E, BetheNoCorr_Linear_MeV_cm(E));
  }
  // --- NIST PSTAR: smooth reference curve + second ratio series ------------
  // The curve is only drawn over the overlap between the plotted range and
  // PSTAR's own tabulated domain — log-log interpolation is not extended
  // past NIST's own data (see InterpolatePstarLogLog).
  TGraph* gPstar = nullptr;
  TGraphErrors* gRatioPstar = nullptr;
  if (havePstar) {
    const double pstarLo = std::max(Elo, pstarData.front().E_MeV);
    const double pstarHi = std::min(Ehi, pstarData.back().E_MeV);
    if (pstarLo < pstarHi) {
      std::vector<double> ex, ey;
      const int nP = 400;
      for (int i = 0; i < nP; ++i) {
        const double E =
            pstarLo * TMath::Power(pstarHi / pstarLo, (double)i / (nP - 1));
        double S = 0.;
        if (InterpolatePstarLogLog(pstarData, E, S)) {
          ex.push_back(E);
          ey.push_back(S * cfg::kDensity);  // mass -> linear stopping power
        }
      }
      if (!ex.empty())
        gPstar = new TGraph((int)ex.size(), ex.data(), ey.data());
    }
    gRatioPstar = new TGraphErrors(n);
    int nRatioP = 0;
    for (int i = 0; i < n; ++i) {
      if (results[i].nExit > 0 && pstarOk[i]) {
        const double pstarLinear = pstarS_MeVcm2_g[i] * cfg::kDensity;
        gRatioPstar->SetPoint(nRatioP, results[i].E,
                              results[i].dedxTotal / pstarLinear);
        gRatioPstar->SetPointError(nRatioP, 0.,
                                   results[i].dedxTotalErr / pstarLinear);
        ++nRatioP;
      }
    }
    gRatioPstar->Set(nRatioP);
  }
  // --- Teoria vs NIST: smooth ratio curve, NO simulation involved ----------
  // Bethe/PSTAR at every energy across the PSTAR-covered part of the plotted
  // range — deterministic (no error bars), so a plain TGraph line rather
  // than TGraphErrors markers, consistent with how gBB/gPstar themselves
  // are drawn as smooth curves rather than per-point markers.
  TGraph* gRatioBethePstar = nullptr;
  if (havePstar) {
    const double pstarLo = std::max(Elo, pstarData.front().E_MeV);
    const double pstarHi = std::min(Ehi, pstarData.back().E_MeV);
    if (pstarLo < pstarHi) {
      std::vector<double> rx, ry;
      const int nR = 400;
      for (int i = 0; i < nR; ++i) {
        const double E =
            pstarLo * TMath::Power(pstarHi / pstarLo, (double)i / (nR - 1));
        double Spstar = 0.;
        if (InterpolatePstarLogLog(pstarData, E, Spstar)) {
          rx.push_back(E);
          ry.push_back(BetheNoCorr_Mass_MeVcm2_g(E) / Spstar);
        }
      }
      if (!rx.empty())
        gRatioBethePstar = new TGraph((int)rx.size(), rx.data(), ry.data());
    }
  }
  // --- Canvas: main panel + ratio panel ------------------------------------
  TCanvas* c = new TCanvas("cScan", "dE/dx vs E", 900, 800);
  c->Divide(1, 2);
  TPad* p1 = (TPad*)c->cd(1);
  p1->SetPad(0., 0.32, 1., 1.);
  p1->SetLogx();
  p1->SetLogy();
  p1->SetBottomMargin(0.02);
  TPad* p2 = (TPad*)c->cd(2);
  p2->SetPad(0., 0., 1., 0.32);
  p2->SetLogx();
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.30);
  p1->cd();
  gBB->SetLineColor(kGray + 2);
  gBB->SetLineWidth(2);
  gBB->SetTitle(Form("Proton stopping power in %s (cut: %s);;-dE/dx (MeV/cm)",
                     cfg::kMaterialName, cfg::kCutTag));
  gBB->GetXaxis()->SetLabelSize(0.);
  gBB->Draw("AL");
  if (gPstar) {
    gPstar->SetLineColor(kGreen + 2);
    gPstar->SetLineStyle(2);
    gPstar->SetLineWidth(2);
    gPstar->Draw("L same");
  }
  gTotal->SetMarkerStyle(20);
  gTotal->SetMarkerColor(kAzure + 2);
  gTotal->SetLineColor(kAzure + 2);
  gTotal->Draw("P same");
  gRestr->SetMarkerStyle(24);
  gRestr->SetMarkerColor(kOrange + 7);
  gRestr->SetLineColor(kOrange + 7);
  gRestr->Draw("P same");
  // Top-right corner: above and to the right of a monotonically decreasing
  // curve on log-log axes there is no data, so the legend cannot overlap
  // points or the curve.
  TLegend* leg = new TLegend(0.44, 0.66, 0.89, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  leg->AddEntry(gBB, "Relativistic Bethe (no corrections)", "l");
  if (gPstar) leg->AddEntry(gPstar, "NIST PSTAR (total, tabulated)", "l");
  leg->AddEntry(gTotal, "Simulation: (E_{in}-E_{out})/track  (~unrestricted)", "p");
  leg->AddEntry(gRestr, "Simulation: E_{dep,primary}/track  (restricted)", "p");
  leg->Draw();
  p2->cd();
  gRatio->SetMarkerStyle(20);
  gRatio->SetMarkerColor(kAzure + 2);
  gRatio->SetLineColor(kAzure + 2);
  gRatio->SetTitle(";proton kinetic energy (MeV);ratio to reference");
  gRatio->GetXaxis()->SetTitleSize(0.10);
  gRatio->GetXaxis()->SetLabelSize(0.09);
  gRatio->GetYaxis()->SetTitleSize(0.09);
  gRatio->GetYaxis()->SetLabelSize(0.08);
  gRatio->GetYaxis()->SetTitleOffset(0.5);
  // The simulation is expected to sit ABOVE the uncorrected Bethe curve
  // (missing shell/Barkas/density corrections), increasingly at low energy;
  // the PSTAR ratios (Geant4/PSTAR and Bethe/PSTAR both include/reflect
  // those corrections) should sit much closer to 1 across the whole band —
  // kept on the SAME 0.9-1.15 range so all series are visually comparable
  // rather than each on its own scale.
  gRatio->GetYaxis()->SetRangeUser(0.9, 1.15);
  gRatio->Draw("AP");
  if (gRatioPstar && gRatioPstar->GetN() > 0) {
    gRatioPstar->SetMarkerStyle(21);
    gRatioPstar->SetMarkerColor(kGreen + 2);
    gRatioPstar->SetLineColor(kGreen + 2);
    gRatioPstar->Draw("P same");
  }
  if (gRatioBethePstar) {
    // Teoria vs NIST: drawn as a smooth GRAY DASHED LINE (matches gBB's
    // color in the main panel) rather than markers, since it's a
    // deterministic curve with no simulation/statistical error attached —
    // visually distinct from the two marker-based sim comparisons above.
    gRatioBethePstar->SetLineColor(kGray + 2);
    gRatioBethePstar->SetLineStyle(2);
    gRatioBethePstar->SetLineWidth(2);
    gRatioBethePstar->Draw("L same");
  }
  TGraph* gUnity = new TGraph(2);
  gUnity->SetPoint(0, Elo, 1.0);
  gUnity->SetPoint(1, Ehi, 1.0);
  gUnity->SetLineColor(kGray + 1);
  gUnity->SetLineStyle(2);
  gUnity->Draw("L same");
  if ((gRatioPstar && gRatioPstar->GetN() > 0) || gRatioBethePstar) {
    // Compact legend in the ratio pad only when there is more than the one
    // (always-present) Geant4-vs-Teoria series to distinguish.
    TLegend* legRatio = new TLegend(0.12, 0.76, 0.68, 0.95);
    legRatio->SetBorderSize(0);
    legRatio->SetFillStyle(0);
    legRatio->SetTextSize(0.070);
    legRatio->SetNColumns(3);
    legRatio->AddEntry(gRatio, "Geant4/Teoria", "p");
    if (gRatioPstar && gRatioPstar->GetN() > 0)
      legRatio->AddEntry(gRatioPstar, "Geant4/NIST", "p");
    if (gRatioBethePstar)
      legRatio->AddEntry(gRatioBethePstar, "Teoria/NIST", "l");
    legRatio->Draw();
  }
  c->SaveAs((outDir + "/dedx_vs_energy.png").c_str());
  c->SaveAs((outDir + "/dedx_vs_energy.pdf").c_str());

  // --- 3D straggling surface ------------------------------------------------
  // Generous margins so the palette, the z-axis title and the axis labels
  // of the LEGO view all fit inside the canvas without clipping or
  // overlapping each other.
  // Up to 6 full digits on axes before switching to scientific notation:
  // suppresses the "x10^3" exponent block that was drawn on top of the
  // event-count scale in both straggling views.
  TGaxis::SetMaxDigits(6);
  TCanvas* cs = new TCanvas("cStrag", "Straggling surface", 1100, 800);
  cs->SetRightMargin(0.18);   // room for the color palette + "events" title
  cs->SetLeftMargin(0.12);
  cs->SetBottomMargin(0.14);  // room for the x-axis title below the labels
  cs->SetTopMargin(0.08);
  hStrag->SetTitleOffset(2.2, "x");
  hStrag->SetTitleOffset(2.0, "y");
  hStrag->SetTitleOffset(1.3, "z");
  hStrag->Draw("LEGO2 Z");
  cs->SaveAs((outDir + "/straggling_3d.png").c_str());
  // Companion 2D color map: same information, easier to read the width
  // evolution quantitatively than the LEGO view.
  TCanvas* cm = new TCanvas("cStragMap", "Straggling map", 1000, 650);
  cm->SetRightMargin(0.15);   // palette
  cm->SetBottomMargin(0.14);  // x-axis title was clipped with the default
  hStrag->GetXaxis()->SetTitleOffset(1.2);
  hStrag->GetYaxis()->SetTitleOffset(1.1);
  hStrag->Draw("COLZ");
  cm->SaveAs((outDir + "/straggling_map.png").c_str());

  // --- Summary CSV for further processing -----------------------------------
  // Linear stopping power in MeV/cm; mass stopping power S in MeV cm2/g
  // (S = linear / density).
  //
  // TOTAL-ESTIMATOR COLUMNS: with the thin-slab thickness rule
  // (t(E) ~ 5% x range(E), applied by generate_energy_scan.py) the total
  // estimator is defined at every sweep energy, so these columns should be
  // complete across the validation band. An empty entry (n_exit == 0) is a
  // DATA-QUALITY FLAG — the run was produced with a thickness incompatible
  // with its energy (e.g. old fixed-5mm data below ~20.5 MeV, or E below
  // the ~0.3 MeV floor where even 1 um exceeds the range) — and such runs
  // must be re-simulated, not compared. Empty fields read as NaN in pandas.
  //
  // rel_err_bethe_vs_total_pct = 100 * (S_total - S_bethe) / S_bethe: the
  // percent deviation of the simulation's ~unrestricted estimator from the
  // uncorrected relativistic Bethe value — the correct pairing, since both
  // are the same physical quantity (the restricted estimator is not: it
  // depends on the production cut). Positive at low energy (missing
  // shell/Barkas corrections), ~0 at 100-1000 MeV, negative at high energy
  // (missing density effect). Empty where the total estimator is undefined.
  //
  // UNCERTAINTIES: bethe_no_corr_err_MeV_cm / S_bethe_no_corr_err_MeV_cm2_g
  // are the ANALYTIC propagated uncertainty on the Bethe reference itself
  // (BetheNoCorr_*_Uncertainty_*, same formula/sources as
  // analytic_solution.ipynb: K, Z/A, m_e, I, +rho for the linear form — NOT
  // a statistical error, since the analytic curve has no "events").
  // rel_err_bethe_vs_total_pct_err combines that analytic uncertainty with
  // the simulation's own statistical error (dedx_total_err_MeV_cm) via the
  // standard uncorrelated propagation formula (RelErrPctUncertainty above).
  // Both are empty where the total estimator is undefined (n_exit == 0).
  //
  // PSTAR COLUMNS: pstar_S_total_MeV_cm2_g / pstar_dedx_total_MeV_cm are the
  // NIST PSTAR "total" (electronic + nuclear) stopping power, log-log
  // interpolated to this run's energy (InterpolatePstarLogLog). Unlike the
  // Bethe columns, PSTAR carries no published uncertainty, so
  // rel_err_pstar_vs_total_pct_err propagates ONLY the simulation's own
  // statistical error (RelErrPctUncertaintyStatOnly) — treat S_pstar as
  // exact, not as having the same error character as bethe_no_corr's.
  // Empty when the total estimator is undefined (n_exit == 0) OR PSTAR data
  // is unavailable/out of range for this energy/material.
  //
  // THEORY VS NIST: bethe_vs_pstar_pct / bethe_vs_pstar_pct_err compare the
  // analytic Bethe curve DIRECTLY against NIST PSTAR (total), with NO
  // dependence on the simulation at all — unlike every other rel_err_*
  // column above, these are defined whenever PSTAR has data at this energy,
  // even for rows where n_exit == 0. Reference (denominator) is PSTAR, so
  // rel_err = 100*(S_bethe - S_pstar)/S_pstar; the propagated uncertainty
  // again uses RelErrPctUncertaintyStatOnly, this time correctly (the
  // "stat-only" side is bethe_no_corr's own analytic uncertainty, S_pstar
  // still treated as exact since NIST publishes none).
  const std::string dedxSummaryPath = outDir + "/dedx_summary.csv";
  std::ofstream out(dedxSummaryPath);
  out << "material,E_MeV,n_events,n_exit,"
         "dedx_restricted_MeV_cm,dedx_restricted_err_MeV_cm,"
         "S_restricted_MeV_cm2_g,"
         "dedx_total_MeV_cm,dedx_total_err_MeV_cm,S_total_MeV_cm2_g,"
         "bethe_no_corr_MeV_cm,bethe_no_corr_err_MeV_cm,"
         "S_bethe_no_corr_MeV_cm2_g,S_bethe_no_corr_err_MeV_cm2_g,"
         "rel_err_bethe_vs_total_pct,rel_err_bethe_vs_total_pct_err,"
         "pstar_S_total_MeV_cm2_g,pstar_dedx_total_MeV_cm,"
         "rel_err_pstar_vs_total_pct,rel_err_pstar_vs_total_pct_err,"
         "rel_err_bethe_vs_pstar_pct,rel_err_bethe_vs_pstar_pct_err,"
         "straggling_rms_MeV,balance_residual_MeV\n";
  // Accumulators for the final cross-material error-summary table (see
  // ErrorSummaryRow / PrintErrorSummaryTable below): mean ABSOLUTE percent
  // deviation for each of the three comparisons, each with its own valid
  // count (the three are independently gated: sim-vs-Bethe needs
  // n_exit > 0; sim-vs-PSTAR needs n_exit > 0 AND PSTAR coverage;
  // Bethe-vs-PSTAR only needs PSTAR coverage).
  double sumAbsBetheVsTotal = 0.; long nAbsBetheVsTotal = 0;
  double sumAbsPstarVsTotal = 0.; long nAbsPstarVsTotal = 0;
  double sumAbsBetheVsPstar = 0.; long nAbsBetheVsPstar = 0;
  for (size_t i = 0; i < results.size(); ++i) {
    const RunResult& r = results[i];
    const double betheLin = BetheNoCorr_Linear_MeV_cm(r.E);
    const double betheLinErr = BetheNoCorr_Linear_Uncertainty_MeV_cm(r.E);
    const double betheMass = BetheNoCorr_Mass_MeVcm2_g(r.E);
    const double betheMassErr = BetheNoCorr_Mass_Uncertainty_MeVcm2_g(r.E);
    out << cfg::kMaterialName << ',' << r.E << ',' << r.n << ',' << r.nExit << ','
        << r.dedxRestricted << ',' << r.dedxRestrictedErr << ','
        << r.dedxRestricted / cfg::kDensity << ',';
    if (r.nExit > 0) {
      const double relErrPct = 100. * (r.dedxTotal - betheLin) / betheLin;
      const double relErrPctErr = RelErrPctUncertainty(
          r.dedxTotal, r.dedxTotalErr, betheLin, betheLinErr);
      out << r.dedxTotal << ',' << r.dedxTotalErr << ','
          << r.dedxTotal / cfg::kDensity << ','
          << betheLin << ',' << betheLinErr << ','
          << betheMass << ',' << betheMassErr << ','
          << relErrPct << ',' << relErrPctErr << ',';
      sumAbsBetheVsTotal += std::fabs(relErrPct);
      ++nAbsBetheVsTotal;
    } else {
      // Primary never exits: total estimator and every comparison against
      // it (Bethe's rel_err columns AND the sim-vs-PSTAR columns below) are
      // undefined (the analytic Bethe value and its own uncertainty are
      // still reported — they don't depend on the simulation). FIVE empty
      // fields here: dedx_total, dedx_total_err, S_total,
      // rel_err_bethe_vs_total_pct, rel_err_bethe_vs_total_pct_err.
      out << ",,," << betheLin << ',' << betheLinErr << ','
          << betheMass << ',' << betheMassErr << ",,,";
    }
    if (r.nExit > 0 && pstarOk[i]) {
      const double pstarMass = pstarS_MeVcm2_g[i];
      const double pstarLin = pstarMass * cfg::kDensity;
      const double relErrPstarPct = 100. * (r.dedxTotal - pstarLin) / pstarLin;
      const double relErrPstarPctErr =
          RelErrPctUncertaintyStatOnly(r.dedxTotal, r.dedxTotalErr, pstarLin);
      out << pstarMass << ',' << pstarLin << ','
          << relErrPstarPct << ',' << relErrPstarPctErr << ',';
      sumAbsPstarVsTotal += std::fabs(relErrPstarPct);
      ++nAbsPstarVsTotal;
    } else {
      out << ",,,,";
    }
    if (pstarOk[i]) {
      // Theory vs NIST: no simulation involved, so this is available
      // whenever PSTAR has data here, regardless of n_exit.
      const double pstarMass = pstarS_MeVcm2_g[i];
      const double relErrBethePstarPct = 100. * (betheMass - pstarMass) / pstarMass;
      const double relErrBethePstarPctErr =
          RelErrPctUncertaintyStatOnly(betheMass, betheMassErr, pstarMass);
      out << relErrBethePstarPct << ',' << relErrBethePstarPctErr << ',';
      sumAbsBetheVsPstar += std::fabs(relErrBethePstarPct);
      ++nAbsBetheVsPstar;
    } else {
      out << ",,";
    }
    out << r.stragglingRms << ',' << r.balanceResidual << '\n';
  }
  out.close();
  printf("\nWritten: %s/dedx_vs_energy.png/.pdf, %s/straggling_3d.png, "
         "%s/straggling_map.png and %s (%d of %zu grid points had data)\n",
         outDir.c_str(), outDir.c_str(), outDir.c_str(), dedxSummaryPath.c_str(),
         n, cfg::kEnergiesMeV.size());

  // --- Record this run in the cross-material error-summary table ----------
  ErrorSummaryRow row;
  row.material = cfg::kMaterialName;
  row.eMinMeV = results.front().E;   // ACTUAL range analyzed: after
  row.eMaxMeV = results.back().E;    // kEMinMeV/kEMaxMeV and available data,
  row.nPoints = (long)results.size();  // not the raw kEnergiesMeV grid.
  row.meanAbsGeant4VsTeoriaPct =
      (nAbsBetheVsTotal > 0) ? sumAbsBetheVsTotal / nAbsBetheVsTotal : 0.;
  row.nGeant4VsTeoria = nAbsBetheVsTotal;
  row.meanAbsGeant4VsNistPct =
      (nAbsPstarVsTotal > 0) ? sumAbsPstarVsTotal / nAbsPstarVsTotal : 0.;
  row.nGeant4VsNist = nAbsPstarVsTotal;
  row.meanAbsTeoriaVsNistPct =
      (nAbsBetheVsPstar > 0) ? sumAbsBetheVsPstar / nAbsBetheVsPstar : 0.;
  row.nTeoriaVsNist = nAbsBetheVsPstar;
  g_errorSummary.push_back(row);
  // Reprints the FULL accumulated table (all materials analyzed so far in
  // this session), not just this run — see PrintErrorSummaryTable() above.
  // Written into THIS run's output directory (see the OUTPUT DIRECTORY note
  // in the header comment) so it always sits next to the most recent data.
  PrintErrorSummaryTable(outDir);
}
// ============================================================================
// Entry point:  root -l -b -q analyze_dedx.C
// ============================================================================
void analyze_dedx(const char* dataDir = ".",
                  const char* pstarDir = cfg::kPstarDir)
{
  printf("=== slab Bethe validation analysis (relativistic Bethe, "
         "no corrections, + NIST PSTAR) ===\n");
  printf("cut tag = %s, thin-slab rule (t <= %g mm), %g <= E <= %g MeV, "
         "grid points = %zu\n", cfg::kCutTag, cfg::kThicknessMaxMM,
         cfg::kEMinMeV, cfg::kEMaxMeV, cfg::kEnergiesMeV.size());
  printf("(material is auto-detected from the data files in '%s'; "
         "NIST PSTAR data read from '%s' — see below)\n\n", dataDir, pstarDir);
  // AnalyzeScan() detects the material (water/aluminium/copper/lead) from
  // the file names present in dataDir, configures cfg::kZoverA/kI_eV/
  // kDensity accordingly, loads the matching NIST PSTAR table from
  // pstarDir (degrading gracefully if not found), and prints the detected
  // material + its analytic Bethe sanity check before doing anything else.
  AnalyzeScan(dataDir, pstarDir);
}