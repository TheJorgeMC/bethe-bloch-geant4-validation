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
// --- Straggling surface (AnalyzeScan) ---------------------------------------
// The deposit of each event is normalized to its run's mean, so runs whose
// absolute deposits differ by orders of magnitude across the sweep share
// one common, comparable y axis. 1.0 = the mean; the shape around it IS
// the straggling distribution.
const int kStragNBinsY = 60;
const double kStragYMin = 0.4;   // Edep / <Edep>
const double kStragYMax = 1.8;
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
// Full sweep: dE/dx vs E against the uncorrected relativistic Bethe curve,
// plus the 3D straggling surface across the sweep
// ============================================================================
void AnalyzeScan(const char* dataDir = ".")
{
  gStyle->SetOptStat(0);
  if (cfg::kEnergiesMeV.empty()) {
    printf("Energy grid is empty — check the path of energy_grid.csv "
           "(read relative to the current directory)\n");
    return;
  }
  // --- Pre-scan: which sweep points actually have data files? -------------
  // Needed to give the straggling TH2 one x bin per available energy.
  std::vector<double> present;
  for (double E : cfg::kEnergiesMeV) {
    if (E < cfg::kEMinMeV) continue;
    TString f = TString::Format(cfg::kNtupleFilePattern, dataDir,
                                EnergyTag(E).Data(), cfg::kCutTag);
    std::ifstream test(f.Data());
    if (test) present.push_back(E);
  }
  if (present.empty()) {
    printf("No data files found in '%s' — check cfg::kNtupleFilePattern and "
           "cfg::kCutTag ('%s')\n", dataDir, cfg::kCutTag);
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
                                EnergyTag(E).Data(), cfg::kCutTag);
    RunData d = LoadNtuple(f.Data());
    RunResult r = AnalyzeRunData(d, E);
    ++ix;
    if (!r.ok) continue;
    results.push_back(r);
    printf("E = %8g MeV : n = %6ld  dE/dx(restr) = %10.4f  "
           "dE/dx(total) = %10.4f  Bethe(no corr) = %10.4f MeV/cm\n",
           E, r.n, r.dedxRestricted, r.dedxTotal,
           BetheNoCorr_Linear_MeV_cm(E));
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
  // Bottom-left corner: empty for a monotonically decreasing curve on
  // log-log axes, so the legend never covers data points.
  TLegend* leg = new TLegend(0.14, 0.06, 0.60, 0.28);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
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
  cs->SaveAs("straggling_3d.png");
  // Companion 2D color map: same information, easier to read the width
  // evolution quantitatively than the LEGO view.
  TCanvas* cm = new TCanvas("cStragMap", "Straggling map", 1000, 650);
  cm->SetRightMargin(0.15);   // palette
  cm->SetBottomMargin(0.14);  // x-axis title was clipped with the default
  hStrag->GetXaxis()->SetTitleOffset(1.2);
  hStrag->GetYaxis()->SetTitleOffset(1.1);
  hStrag->Draw("COLZ");
  cm->SaveAs("straggling_map.png");

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
  printf("\nWritten: dedx_vs_energy.png/.pdf, straggling_3d.png, "
         "straggling_map.png and dedx_summary.csv "
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