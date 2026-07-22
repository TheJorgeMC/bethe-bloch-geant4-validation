// ============================================================================
// Run.hh
// G4Run subclass that accumulates, thread-safely (each worker thread owns its
// local instance; G4 calls Merge() on the master's), the sums needed for the
// end-of-run statistics: mean and standard deviation of the energy deposit
// (straggling), per-channel energy balance, and the comparison against
// G4EmCalculator in the style of the TestEm0 report.
// ============================================================================
#ifndef Run_hh
#define Run_hh 1

#include "G4Run.hh"
#include "globals.hh"

class DetectorConstruction;
class G4ParticleDefinition;

class Run : public G4Run
{
 public:
  explicit Run(const DetectorConstruction* det);
  ~Run() override = default;

  // Registered by RunAction::BeginOfRunAction on the worker threads from the
  // G4ParticleGun (TestEm examples pattern).
  void SetPrimary(const G4ParticleDefinition* particle, G4double energy);

  // Registered by EventAction::EndOfEventAction, once per event.
  //   exitEnergy < 0 means the primary did NOT exit the slab.
  void AddEvent(G4double edepPrimary, G4double edepSecondary,
                G4double escapedSecondary, G4double exitEnergy,
                G4double trackLength);

  void Merge(const G4Run* aRun) override;

  // Console report (called by the master only, in EndOfRunAction).
  void EndOfRun() const;

 private:
  const DetectorConstruction* fDetector;

  const G4ParticleDefinition* fParticle = nullptr;
  G4double fEkin = 0.;

  G4double fSumEdepTotal = 0.;   // primary + secondaries (per event)
  G4double fSumEdepTotal2 = 0.;  // sum of squares, for the straggling
  G4double fSumEdepPrimary = 0.;
  G4double fSumEdepSecondary = 0.;
  G4double fSumEscapedSecondary = 0.;
  G4double fSumExitEnergy = 0.;  // only events where the primary exits
  G4long fNExit = 0;             // number of events with an exiting primary
  G4double fSumTrackLength = 0.;
  G4long fNEvents = 0;
};

#endif
