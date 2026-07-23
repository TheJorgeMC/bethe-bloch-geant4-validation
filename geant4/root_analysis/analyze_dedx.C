// ============================================================================
// analyze_dedx.C — ROOT analysis macro for the slab Bethe-Bloch validation
//
// Reads the merged CSV output of the Geant4 application (one
// "<base>_nt_slab.csv" ntuple per run and one "<base>_h1_DepthEdep.csv"
// depth-dose histogram) and produces:
//
//   1) Per-run statistics: mean/rms energy deposit (straggling), restricted
//      and total dE/dx estimators, per-event energy-balance residual.
//   2) The dE/dx vs energy curve of the full sweep, compared against the
//      PURE RELATIVISTIC BETHE formula (NO corrections) evaluated in this
//      macro — see the note below.
//   3) The Bragg / depth-dose curve from the DepthEdep histogram.
//
// ANALYTIC REFERENCE (consistent with analytic_solution.ipynb and the
// "Relativistic Bethe" section of Bethe_Full_Derivation.pdf):
//
//   S(T) = K z^2 (Z/A) (1/beta^2) [ ln(2 m_e c^2 beta^2 gamma^2 / I) - beta^2 ]
//
// deliberately EXCLUDING shell (C/Z), Barkas (z^3), Bloch and density
// (Sternheimer delta) corrections, and using the heavy-projectile
// approximation Tmax ~ 2 m_e c^2 beta^2 gamma^2 inside the logarithm.
// Constants are CODATA 2022 / PDG; water Z/A follows the Bragg additivity
// rule with NIST STAR mass fractions; I = 78(2) eV from ICRU 90.
// EXPECTED BEHAVIOR: this curve systematically UNDERESTIMATES both NIST
// PSTAR and the Geant4 simulation, increasingly below ~10 MeV — that is a
// documented consequence of the project scope, NOT a bug.
//
// FILE-NAME TAGS: energies (and the cut) are 'p'-encoded in file names,
// because G4AnalysisManager treats everything after the last '.' as a file
// extension: 0.001 MeV -> "dedx_0p001MeV_cut0p01mm_nt_slab.csv". The
// EnergyTag() helper below must stay consistent with etag() in
// generate_energy_scan.py.
//
// Usage (from the directory containing the CSV files):
//   root -l -b -q analyze_dedx.C                      # full sweep + plots
//   root -l 'analyze_dedx.C("path/to/data")'          # data in another dir
//   root -l
//     .L analyze_dedx.C
//     AnalyzeRun("dedx_150MeV_cut0p01mm_nt_slab.csv", 150.0);  // single run
//     PlotBragg("dedx_150MeV_cut0p01mm_h1_DepthEdep.csv");
//
// The configuration block below (energy grid CSV, cut tag, slab thickness,
// material constants) must match the macros used to produce the data.
// ============================================================================
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"
#include "TPad.h"
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
namespace cfg
{
// Minimum energy included in the analysis (MeV). Points below this are
// skipped: the uncorrected Bethe formula is outside its Born-approximation
// validity range there (see the Sommerfeld-parameter discussion in the
// derivation document) and sub-range protons stop inside the slab anyway.
const double kEMinMeV = 3.0;
// Sweep energies in MeV: NIST PSTAR preferred-number grid, same file used by
// generate_energy_scan.py (path relative to where root is executed).
const std::vector<double> kEnergiesMeV =
    readFirstColumnFromCSV("../../nist_data/energy_grid.csv");
// Cut tag in the data file names — must match CUT_TAG in
// generate_energy_scan.py ("cut0p01mm" for the 0.01 mm fine cut;
// "cut1mm" for the 1 mm default).
const char* kCutTag = "cut0p01mm";
// Ntuple file pattern: dedx_<Etag>MeV_<cutTag>_nt_slab.csv
// (%s placeholders: dataDir, EnergyTag(E), kCutTag)
const char* kNtupleFilePattern = "%s/dedx_%sMeV_%s_nt_slab.csv";
// Slab thickness in mm (must match /absorber/thickness in the macros).
// Only used for reference printing; the dE/dx estimators divide by the
// actual per-event track length from the ntuple, which is more accurate.
const double kThicknessMM = 5.0;
// --- Absorber material constants: liquid water -----------------------------
// Kept numerically identical to the WATER Material in analytic_solution.ipynb.
// Z/A by the Bragg additivity rule (ICRU 37/49, Sec. 2.5.2 Eq. 2.22) with
// NIST STAR mass fractions (matno 276) and IUPAC 2021 atomic weights:
//   H: w = 0.111894, A = 1.0080(2);  O: w = 0.888106, A = 15.999(1).
const double kWaterWH = 0.111894, kWaterAH = 1.0080, kWaterZH = 1.0;
const double kWaterWO = 0.888106, kWaterAO = 15.999, kWaterZO = 8.0;
const double kZoverA =
    kWaterWH * (kWaterZH / kWaterAH) + kWaterWO * (kWaterZO / kWaterAO);
const double kDensity = 1.0;      // g/cm3 (PubChem; as in the notebook)
const double kI_eV = 78.0;        // mean excitation energy, ICRU 90: 78(2) eV
// --- Physical constants (CODATA 2022 / PDG, as in the notebook) -------------
const double kMp = 938.27208943;    // proton rest mass energy, MeV
const double kMe = 0.51099895069;   // electron rest mass energy, MeV
const double kK = 0.307075;         // 4 pi N_A r_e^2 m_e c^2, MeV cm2 / mol
const double kEvToMeV = 1e-6;       // eV -> MeV
}  // namespace cfg
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
struct RunResult
{
  double E = 0.;                    // nominal beam energy (MeV)
  // Linear stopping power estimators in MeV/cm (ntuple lengths are in mm;
  // the mm -> cm conversion is applied in AnalyzeRunData). Mass stopping
  // power (MeV cm2/g) = linear / density.
  double dedxRestricted = 0.;       // <Edep_primary / track_length>  (MeV/cm)
  double dedxRestrictedErr = 0.;
  double dedxTotal = 0.;            // <(E_inc - E_exit) / track_length> (MeV/cm)
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
  RunData d = LoadNtuple(csvFile);
  if (!d.ok) return;
  RunResult r = AnalyzeRunData(d, Enominal);
  printf("\n--- %s ---\n", csvFile);
  printf("  events                  : %ld\n", r.n);
  printf("  exiting primaries       : %ld / %ld\n", r.nExit, r.n);
  printf("  dE/dx restricted (sim)  : %.4f +- %.4f MeV/cm  "
         "(S = %.4f MeV cm2/g)\n",
         r.dedxRestricted, r.dedxRestrictedErr,
         r.dedxRestricted / cfg::kDensity);
  printf("  dE/dx total      (sim)  : %.4f +- %.4f MeV/cm  "
         "(S = %.4f MeV cm2/g)\n",
         r.dedxTotal, r.dedxTotalErr, r.dedxTotal / cfg::kDensity);
  printf("  Bethe, no corrections   : %.4f MeV/cm  (S = %.4f MeV cm2/g)\n",
         BetheNoCorr_Linear_MeV_cm(Enominal),
         BetheNoCorr_Mass_MeVcm2_g(Enominal));
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
  // Straggling-shape fits.
  // The energy-loss distribution in a slab interpolates between the two
  // classical limits of Landau-Vavilov theory, governed by the kappa
  // parameter (mean energy loss / Tmax): thin absorber (kappa << 1) ->
  // asymmetric Landau with a delta-ray tail; thick absorber (kappa >> 1)
  // -> Gaussian (Bohr straggling). Fitting BOTH and comparing chi2/ndf
  // shows which regime the run sits in. ROOT's built-in "landau" TF1 is
  // the standard proxy for the Landau-Vavilov shape (MPV and width scale
  // parameter; exact Vavilov fitting is rarely warranted at this level).
  // ---------------------------------------------------------------------
  TF1* fGaus = new TF1(Form("fGaus_%s", tag.Data()), "gaus", lo, hi);
  fGaus->SetParameters(h->GetMaximum(), mean, rms);
  h->Fit(fGaus, "Q0R");  // Q: quiet, 0: don't auto-draw, R: fit range

  TF1* fLandau = new TF1(Form("fLandau_%s", tag.Data()), "landau", lo, hi);
  // Seeds: MPV slightly below the mean (Landau mean > MPV), width ~ rms/4.
  fLandau->SetParameters(h->GetMaximum(), mean - 0.2 * rms, 0.25 * rms);
  h->Fit(fLandau, "Q0R+");  // +: add to the fit list, keep the Gaussian

  const double gausChi2 =
      (fGaus->GetNDF() > 0) ? fGaus->GetChisquare() / fGaus->GetNDF() : 0.;
  const double landauChi2 =
      (fLandau->GetNDF() > 0) ? fLandau->GetChisquare() / fLandau->GetNDF() : 0.;
  printf("  Gaussian fit            : mean = %.4f +- %.4f MeV, "
         "sigma = %.4f +- %.4f MeV, chi2/ndf = %.2f\n",
         fGaus->GetParameter(1), fGaus->GetParError(1),
         fGaus->GetParameter(2), fGaus->GetParError(2), gausChi2);
  printf("  Landau-Vavilov fit      : MPV  = %.4f +- %.4f MeV, "
         "width = %.4f +- %.4f MeV, chi2/ndf = %.2f\n",
         fLandau->GetParameter(1), fLandau->GetParError(1),
         fLandau->GetParameter(2), fLandau->GetParError(2), landauChi2);

  TCanvas* c = new TCanvas(Form("cEdep_%s", tag.Data()), "Energy deposit", 800, 600);
  h->SetLineColor(kAzure + 2);
  h->SetLineWidth(2);
  h->SetStats(false);
  h->Draw("hist");
  fGaus->SetLineColor(kOrange + 7);
  fGaus->SetLineWidth(2);
  fGaus->Draw("same");
  fLandau->SetLineColor(kGreen + 2);
  fLandau->SetLineWidth(2);
  fLandau->SetLineStyle(2);
  fLandau->Draw("same");

  TLegend* legFit = new TLegend(0.14, 0.62, 0.52, 0.88);
  legFit->SetBorderSize(0);
  legFit->SetFillStyle(0);
  legFit->AddEntry(h, "Simulation", "l");
  legFit->AddEntry(fGaus,
                   Form("Gaussian: #mu=%.3f, #sigma=%.3f (#chi^{2}/ndf=%.1f)",
                        fGaus->GetParameter(1), fGaus->GetParameter(2), gausChi2),
                   "l");
  legFit->AddEntry(fLandau,
                   Form("Landau: MPV=%.3f, w=%.3f (#chi^{2}/ndf=%.1f)",
                        fLandau->GetParameter(1), fLandau->GetParameter(2),
                        landauChi2),
                   "l");
  legFit->Draw();

  // 'p'-encoded tag in the output name too, for consistent sorting/globbing.
  c->SaveAs(Form("edep_distribution_%sMeV.png", tag.Data()));
}
// ============================================================================
// Depth-dose (Bragg) curve from the Geant4 CSV histogram file
// (g4tools wcsv format: '#' header lines with "#axis fixed N min max",
//  a column-title line, then N+2 rows "entries,Sw,Sw2,Sxw,Sx2w" where the
//  first and last rows are the under/overflow bins)
// ============================================================================
void PlotBragg(const char* h1File)
{
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
  c->SaveAs("bragg_curve.png");
  printf("Bragg curve written to bragg_curve.png (%d bins, z in [%g, %g] mm)\n",
         nbins, zmin, zmax);
}
// ============================================================================
// Full sweep: dE/dx vs E against the uncorrected relativistic Bethe curve
// ============================================================================
void AnalyzeScan(const char* dataDir = ".")
{
  gStyle->SetOptStat(0);
  if (cfg::kEnergiesMeV.empty()) {
    printf("Energy grid is empty — check the path of energy_grid.csv "
           "(read relative to the current directory)\n");
    return;
  }
  std::vector<RunResult> results;
  for (double E : cfg::kEnergiesMeV) {
    if (E < cfg::kEMinMeV) continue;
    // 'p'-encoded energy tag in the file name (0.001 -> 0p001), matching
    // etag() in generate_energy_scan.py.
    TString f = TString::Format(cfg::kNtupleFilePattern, dataDir,
                                EnergyTag(E).Data(), cfg::kCutTag);
    RunData d = LoadNtuple(f.Data());
    RunResult r = AnalyzeRunData(d, E);
    if (r.ok) {
      results.push_back(r);
      printf("E = %8g MeV : n = %6ld  dE/dx(restr) = %10.4f  "
             "dE/dx(total) = %10.4f  Bethe(no corr) = %10.4f MeV/cm\n",
             E, r.n, r.dedxRestricted, r.dedxTotal,
             BetheNoCorr_Linear_MeV_cm(E));
    }
  }
  if (results.empty()) {
    printf("No data files found in '%s' — check cfg::kNtupleFilePattern and "
           "cfg::kCutTag ('%s')\n", dataDir, cfg::kCutTag);
    return;
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
    // Ratio only where the total estimator exists (primaries that exit the
    // slab): below ~20 MeV the proton stops inside 5 mm of water and
    // (E_in - E_out) is undefined.
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
  gBB->SetTitle(Form("Proton stopping power in water (cut: %s);;-dE/dx (MeV/cm)", cfg::kCutTag));
  gBB->GetXaxis()->SetLabelSize(0.);
  gBB->Draw("AL");
  gTotal->SetMarkerStyle(20);
  gTotal->SetMarkerColor(kAzure + 2);
  gTotal->SetLineColor(kAzure + 2);
  gTotal->Draw("P same");
  gRestr->SetMarkerStyle(24);
  gRestr->SetMarkerColor(kOrange + 7);
  gRestr->SetLineColor(kOrange + 7);
  gRestr->Draw("P same");
  TLegend* leg = new TLegend(0.42, 0.65, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(gBB, "Relativistic Bethe (no corrections)", "l");
  leg->AddEntry(gTotal, "Simulation: (E_{in}-E_{out})/track  (~unrestricted)", "p");
  leg->AddEntry(gRestr, "Simulation: E_{dep,primary}/track  (restricted)", "p");
  leg->Draw();
  p2->cd();
  gRatio->SetMarkerStyle(20);
  gRatio->SetMarkerColor(kAzure + 2);
  gRatio->SetLineColor(kAzure + 2);
  gRatio->SetTitle(";proton kinetic energy (MeV);sim(total)/Bethe(no corr)");
  gRatio->GetXaxis()->SetTitleSize(0.10);
  gRatio->GetXaxis()->SetLabelSize(0.09);
  gRatio->GetYaxis()->SetTitleSize(0.09);
  gRatio->GetYaxis()->SetLabelSize(0.08);
  gRatio->GetYaxis()->SetTitleOffset(0.5);
  // The simulation is expected to sit ABOVE the uncorrected Bethe curve
  // (missing shell/Barkas/density corrections), increasingly at low energy.
  gRatio->GetYaxis()->SetRangeUser(0.9, 1.15);
  gRatio->Draw("AP");
  TGraph* gUnity = new TGraph(2);
  gUnity->SetPoint(0, Elo, 1.0);
  gUnity->SetPoint(1, Ehi, 1.0);
  gUnity->SetLineColor(kGray + 1);
  gUnity->SetLineStyle(2);
  gUnity->Draw("L same");
  c->SaveAs("dedx_vs_energy.png");
  c->SaveAs("dedx_vs_energy.pdf");
  // --- Summary CSV for further processing -----------------------------------
  // Linear stopping power in MeV/cm; mass stopping power S in MeV cm2/g
  // (S = linear / density).
  std::ofstream out("dedx_summary.csv");
  out << "E_MeV,n_events,n_exit,"
         "dedx_restricted_MeV_cm,dedx_restricted_err_MeV_cm,"
         "S_restricted_MeV_cm2_g,"
         "dedx_total_MeV_cm,dedx_total_err_MeV_cm,S_total_MeV_cm2_g,"
         "bethe_no_corr_MeV_cm,S_bethe_no_corr_MeV_cm2_g,"
         "straggling_rms_MeV,balance_residual_MeV\n";
  for (const RunResult& r : results) {
    out << r.E << ',' << r.n << ',' << r.nExit << ','
        << r.dedxRestricted << ',' << r.dedxRestrictedErr << ','
        << r.dedxRestricted / cfg::kDensity << ','
        << r.dedxTotal << ',' << r.dedxTotalErr << ','
        << r.dedxTotal / cfg::kDensity << ','
        << BetheNoCorr_Linear_MeV_cm(r.E) << ','
        << BetheNoCorr_Mass_MeVcm2_g(r.E) << ','
        << r.stragglingRms << ',' << r.balanceResidual << '\n';
  }
  out.close();
  printf("\nWritten: dedx_vs_energy.png/.pdf and dedx_summary.csv "
         "(%d of %zu grid points had data)\n",
         n, cfg::kEnergiesMeV.size());
}
// ============================================================================
// Entry point:  root -l -b -q analyze_dedx.C
// ============================================================================
void analyze_dedx(const char* dataDir = ".")
{
  printf("=== slab Bethe validation analysis (relativistic Bethe, "
         "no corrections) ===\n");
  printf("cut tag = %s, I = %g eV, Z/A = %.6f, rho = %g g/cm3, "
         "thickness = %g mm, E >= %g MeV, grid points = %zu\n\n",
         cfg::kCutTag, cfg::kI_eV, cfg::kZoverA, cfg::kDensity,
         cfg::kThicknessMM, cfg::kEMinMeV, cfg::kEnergiesMeV.size());
  AnalyzeScan(dataDir);
}