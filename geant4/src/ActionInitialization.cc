// ============================================================================
// ActionInitialization.cc
// ============================================================================
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(const DetectorConstruction* det)
  : G4VUserActionInitialization(), fDetector(det)
{
}

void ActionInitialization::BuildForMaster() const
{
  // The master generates no primaries in MT mode: RunAction receives nullptr
  // and the master's Run obtains particle/energy through Run::Merge()
  // (TestEm pattern).
  SetUserAction(new RunAction(fDetector, nullptr));
}

void ActionInitialization::Build() const
{
  auto* primary = new PrimaryGeneratorAction(fDetector);
  SetUserAction(primary);

  SetUserAction(new RunAction(fDetector, primary));

  auto* eventAction = new EventAction(fDetector);
  SetUserAction(eventAction);

  SetUserAction(new SteppingAction(fDetector, eventAction));
}
