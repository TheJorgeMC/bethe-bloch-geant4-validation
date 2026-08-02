// ============================================================================
// RunAction.cc
// ============================================================================

#include "RunAction.hh"

#include "AnalysisDefs.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

namespace
{
// Default output file name; overridden with /analysis/setFileName.
const G4String kDefaultFileName = "slab_out";

// Default file type. CSV is directly readable with pandas; for a natively
// merged single file use a .root extension in /analysis/setFileName
// (e.g. dedx_150MeV.root) — see README. For CSV in MT mode, the per-thread
// files are merged by MergeCsvNtupleFiles() below.
const G4String kDefaultFileType = "csv";

// Name of the ntuple booked in Book() (used to locate the CSV files).
const G4String kNtupleName = "slab";

// Returns true when the analysis file name asks for ROOT output, which is the
// only format with native cross-thread ntuple merging.
G4bool WantsRootOutput(const std::string& fileName)
{
  const std::string ext = ".root";
  return fileName.size() > ext.size() &&
         fileName.compare(fileName.size() - ext.size(), ext.size(), ext) == 0;
}
}  // namespace

RunAction::RunAction(const DetectorConstruction* det,
                     const PrimaryGeneratorAction* primary)
  : G4UserRunAction(), fDetector(det), fPrimary(primary)
{
  Book();
}

void RunAction::Book()
{
  auto* am = G4AnalysisManager::Instance();
  am->SetDefaultFileType(kDefaultFileType);
  am->SetFileName(kDefaultFileName);
  am->SetVerboseLevel(0);
  // NOTE: SetNtupleMerging() is deliberately NOT called here. Native merging
  // only exists for ROOT output, and requesting it with CSV makes Geant4 emit
  // "Ntuple merging is not available with csv output" once per thread per run.
  // The decision needs the actual file name, which is only known after
  // /analysis/setFileName has been processed, so it is taken in
  // BeginOfRunAction instead. For CSV, MergeCsvNtupleFiles() does the job.

  // --- H1 kH1Depth: depth dose. Created with provisional binning and
  // rebinned in BeginOfRunAction with the current layer count and thickness
  // (SetH1 allows redefining bins between runs). Histograms ARE merged
  // natively across threads for every output format, so this one always
  // yields a single file.
  const G4int h1 = am->CreateH1(
      "DepthEdep", "Energy deposited vs depth;z (mm);Edep (MeV)", 1, 0., 1.);
  if (h1 != analysis::kH1Depth) {
    G4Exception("RunAction::Book()", "SlabAnalysis001", FatalException,
                "The DepthEdep histogram id does not match kH1Depth");
  }

  // --- Per-event ntuple (one row per primary proton) ---
  am->CreateNtuple(kNtupleName, "Per-event energy balance (MeV, mm)");
  am->CreateNtupleDColumn("E_incident");         // analysis::kColEIncident
  am->CreateNtupleDColumn("Edep_primary");       // analysis::kColEdepPrimary
  am->CreateNtupleDColumn("Edep_secondary");     // analysis::kColEdepSecondary
  am->CreateNtupleDColumn("E_escaped_secondary");// analysis::kColEEscapedSec
  am->CreateNtupleDColumn("E_exit_primary");     // analysis::kColEExitPrimary
  am->CreateNtupleDColumn("track_length_primary");// analysis::kColTrackLenPrimary
  am->FinishNtuple();
}

G4Run* RunAction::GenerateRun()
{
  fRun = new Run(fDetector);
  return fRun;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  // On worker threads, register particle and energy from the gun into the
  // local Run; Merge() propagates them to the master's Run (TestEm pattern).
  if (fPrimary != nullptr && fRun != nullptr) {
    fRun->SetPrimary(fPrimary->GetGun()->GetParticleDefinition(),
                     fPrimary->GetGun()->GetParticleEnergy());
  }

  auto* am = G4AnalysisManager::Instance();

  // Rebin the depth histogram with the current geometry: one bin per layer
  // when the slab is segmented; a single bin otherwise.
  const G4int nLayers = fDetector->GetNumberOfLayers();
  const G4double thicknessMM = fDetector->GetThickness() / mm;
  am->SetH1(analysis::kH1Depth, nLayers, 0., thicknessMM);

  const std::string fileName = am->GetFileName();

  // Enable native cross-thread ntuple merging only when the output really is
  // ROOT; with CSV it is ignored and only produces a warning per thread. Must
  // be set before OpenFile().
  am->SetNtupleMerging(WantsRootOutput(fileName));

  // G4AnalysisManager does NOT create directories: the name given to
  // /analysis/setFileName is handed straight to the file-open call, and the
  // path is resolved relative to the CWD of the executable, not to wherever
  // generate_energy_scan.py happened to run. In the multi-material sweep the
  // name carries a "material_data/<key>/" prefix, so the directory must be
  // created here or the open silently fails and every FillNtupleDColumn()
  // then raises "Ntuple 0 does not exist."
  EnsureOutputDirectory(fileName);

  am->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* am = G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();

  // The console report with statistics merged from all threads only makes
  // sense on the master. The master's EndOfRunAction runs after every worker
  // has finished and closed its file, so this is also the safe place to
  // merge the per-thread CSV ntuple files into a single one.
  if (IsMaster()) {
    MergeCsvNtupleFiles(kNtupleName);
    if (fRun != nullptr) fRun->EndOfRun();
  }
}

void RunAction::EnsureOutputDirectory(const G4String& fileName) const
{
  const fs::path outPath = std::string(fileName);
  if (!outPath.has_parent_path()) return;
  const fs::path dir = outPath.parent_path();
  if (dir.empty()) return;

  std::error_code ec;
  fs::create_directories(dir, ec);
  // create_directories() reports no error when the directory already exists,
  // so a non-empty ec here is a genuine failure (permissions, a file with the
  // same name, ...). Every worker thread calls this; the operation is
  // idempotent and the error_code overload never throws on a benign race.
  if (ec && !fs::is_directory(dir)) {
    G4cerr << "### RunAction: could not create output directory "
           << dir.string() << " : " << ec.message() << G4endl;
  }
}

void RunAction::MergeCsvNtupleFiles(const G4String& ntupleName) const
{
  // Base name as set via /analysis/setFileName. If the user requested a
  // non-CSV format through an explicit extension (.root, .hdf5, .xml),
  // Geant4 merges natively (or writes a single file) and there is nothing
  // to do here. A trailing ".csv" extension, if given, is stripped because
  // Geant4 builds the ntuple file names from the bare base name.
  std::string base = G4AnalysisManager::Instance()->GetFileName();
  const auto dot = base.find_last_of('.');
  if (dot != std::string::npos) {
    const std::string ext = base.substr(dot + 1);
    if (ext == "root" || ext == "hdf5" || ext == "xml") return;
    if (ext == "csv") base.erase(dot);
  }

  // The base name may carry a directory prefix ("material_data/<key>/dedx_...").
  // Geant4 writes the per-thread files next to it, so the search has to happen
  // in that directory and the prefix match has to be done against the bare file
  // name — comparing a path-carrying prefix against filename() would never
  // match and the merge would silently do nothing.
  const fs::path basePath(base);
  const fs::path dir =
      basePath.has_parent_path() && !basePath.parent_path().empty()
          ? basePath.parent_path()
          : fs::current_path();
  const std::string stem = basePath.filename().string();

  // Geant4 names the worker CSV files "<stem>_nt_<ntuple>_t<i>.csv".
  const std::string prefix = stem + "_nt_" + std::string(ntupleName) + "_t";
  const std::string suffix = ".csv";

  struct ThreadFile
  {
    long id;
    fs::path path;
  };
  std::vector<ThreadFile> files;

  std::error_code ec;
  for (const auto& entry : fs::directory_iterator(dir, ec)) {
    if (ec) break;
    if (!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();
    if (name.rfind(prefix, 0) != 0) continue;
    if (name.size() <= prefix.size() + suffix.size()) continue;
    if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
      continue;
    const std::string idStr =
        name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
    if (idStr.empty() ||
        !std::all_of(idStr.begin(), idStr.end(),
                     [](unsigned char c) { return std::isdigit(c); }))
      continue;
    files.push_back({std::stol(idStr), entry.path()});
  }

  // Sequential mode (single "<base>_nt_<ntuple>.csv", no _t suffix) or no
  // output at all: nothing to merge.
  if (files.empty()) return;

  std::sort(files.begin(), files.end(),
            [](const ThreadFile& a, const ThreadFile& b) { return a.id < b.id; });

  const fs::path mergedPath =
      dir / (stem + "_nt_" + std::string(ntupleName) + ".csv");
  std::ofstream out(mergedPath, std::ios::trunc);
  if (!out) {
    G4cerr << "### RunAction: could not open " << mergedPath.string()
           << " for writing; per-thread files were left untouched" << G4endl;
    return;
  }

  // Keep the header block (lines starting with '#': column declarations)
  // from the first file only; append the data rows of every file in thread
  // order. Row order across threads is irrelevant: events are statistically
  // independent.
  bool headerWritten = false;
  std::size_t nRows = 0;
  for (const auto& tf : files) {
    std::ifstream in(tf.path);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      if (line[0] == '#') {
        if (!headerWritten) out << line << '\n';
        continue;
      }
      out << line << '\n';
      ++nRows;
    }
    headerWritten = true;
  }

  out.close();

  ec.clear();
  for (const auto& tf : files) fs::remove(tf.path, ec);

  G4cout << "### RunAction: merged " << files.size()
         << " per-thread ntuple files (" << nRows << " rows) into "
         << mergedPath.string() << G4endl;
}