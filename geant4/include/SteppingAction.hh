// ============================================================================
// SteppingAction.hh
// Manual step-by-step scoring inside the slab: separates the energy
// deposited by the primary proton (ParentID == 0) from that of secondaries
// (delta rays and their progeny), accounts for the energy that secondaries
// carry out of the slab, and records the primary's exit energy and track
// length. This is the mechanism that reconstructs the balance
//   dE/dx(unrestricted) = dE/dx(restricted) + energy of explicit deltas
// from the design document (section 5.3).
// ============================================================================
#ifndef SteppingAction_hh
#define SteppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;

class SteppingAction : public G4UserSteppingAction
{
 public:
  SteppingAction(const DetectorConstruction* det, EventAction* eventAction);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  const DetectorConstruction* fDetector;
  EventAction* fEventAction;
};

#endif
