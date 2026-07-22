// ============================================================================
// RunAction.hh
// Opens/closes the G4AnalysisManager file (CSV by default; ROOT if the file
// name carries a .root extension), defines the per-event ntuple and the
// depth-dose histogram, and delegates the end-of-run console report to
// Run::EndOfRun() (master only).
//
// CSV in multi-threaded mode: Geant4 has no native ntuple merging for the
// CSV format (only ROOT/HDF5), so each worker thread writes its own
// "<name>_nt_slab_t<i>.csv" file. To deliver a single output file, the
// master merges all per-thread ntuple files into "<name>_nt_slab.csv" at
// the end of each run and removes the per-thread files (MergeCsvNtupleFiles).
// ============================================================================
#ifndef RunAction_hh
#define RunAction_hh 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;
class Run;

class RunAction : public G4UserRunAction
{
 public:
  // primary == nullptr on the master thread (TestEm pattern: the master has
  // no primary generator in multi-threaded mode).
  RunAction(const DetectorConstruction* det,
            const PrimaryGeneratorAction* primary);
  ~RunAction() override = default;

  G4Run* GenerateRun() override;
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  void Book();  // single definition of ntuple and histograms

  // Merges the per-thread CSV ntuple files ("<base>_nt_<ntuple>_t<i>.csv")
  // into a single "<base>_nt_<ntuple>.csv" and deletes the per-thread files.
  // No-op in sequential mode or for ROOT output.
  void MergeCsvNtupleFiles(const G4String& ntupleName) const;

  const DetectorConstruction* fDetector;
  const PrimaryGeneratorAction* fPrimary;
  Run* fRun = nullptr;
};

#endif
