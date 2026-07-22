// ============================================================================
// SteppingAction.cc
// ============================================================================
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

SteppingAction::SteppingAction(const DetectorConstruction* det,
                               EventAction* eventAction)
  : G4UserSteppingAction(), fDetector(det), fEventAction(eventAction)
{
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Only steps that OCCUR inside the slab matter: the pre-step point's
  // volume defines where the step's energy was deposited.
  const G4StepPoint* prePoint = step->GetPreStepPoint();
  const G4LogicalVolume* preLV =
      prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (!fDetector->IsAbsorber(preLV)) return;

  G4Track* track = step->GetTrack();
  const G4bool isPrimary = (track->GetParentID() == 0);

  // Layer index for depth scoring. With G4PVReplica the replica number is
  // the copy number of the touchable at depth 0. With an unsegmented slab
  // (numberOfLayers == 1) the copy number is 0.
  const G4int layer = (fDetector->GetNumberOfLayers() > 1)
                          ? prePoint->GetTouchableHandle()->GetCopyNumber()
                          : 0;

  // --------------------------------------------------------------------------
  // 1) Energy deposited in this step.
  // PHYSICS: for the primary, GetTotalEnergyDeposit() contains the
  // CONTINUOUS energy loss of the ionisation process — i.e. the RESTRICTED
  // dE/dx with T_cut set by the region's production cut — plus any
  // below-threshold discrete deposits. The energy handed to EXPLICIT delta
  // rays is NOT included here: it travels with the secondary track and will
  // show up in the secondaries counter (if the delta stops inside the slab)
  // or in the escaped-energy counter (if it leaves).
  // --------------------------------------------------------------------------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) {
    if (isPrimary) {
      fEventAction->AddPrimaryEdep(edep, layer);
    } else {
      fEventAction->AddSecondaryEdep(edep, layer);
    }
  }

  // 2) Primary track length inside the slab (to normalize by the actual
  // path instead of the geometric thickness when multiple scattering is
  // appreciable).
  if (isPrimary) {
    fEventAction->AddPrimaryTrackLength(step->GetStepLength());
  }

  // --------------------------------------------------------------------------
  // 3) Boundary crossing out of the slab.
  // If the post-step point sits on a geometric boundary and the next volume
  // no longer belongs to the slab (note: moving from one replicated layer
  // to the next is still slab; IsAbsorber covers that), the track is
  // leaving.
  // --------------------------------------------------------------------------
  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    const G4VPhysicalVolume* nextPV = track->GetNextVolume();
    const G4LogicalVolume* nextLV =
        (nextPV != nullptr) ? nextPV->GetLogicalVolume() : nullptr;
    if (nextLV == nullptr || !fDetector->IsAbsorber(nextLV)) {
      const G4double ekin = track->GetKineticEnergy();
      if (isPrimary) {
        // Primary kinetic energy on leaving the slab (through the exit face
        // or, after heavy scattering, through a side: recorded either way;
        // the ntuple allows filtering by track length).
        fEventAction->SetPrimaryExitEnergy(ekin);
      } else if (ekin > 0.) {
        // Energy that a secondary (typically a delta ray or a
        // bremsstrahlung/fluorescence photon) CARRIES OUT of the slab
        // without depositing it. This is the component to add back to the
        // local deposit to recover the unrestricted dE/dx (balance of
        // section 5.3).
        // NOTE: if the same track re-entered the slab and left again it
        // would be counted twice; with this geometry (slab surrounded by
        // air) the effect is negligible and accepted as an approximation.
        fEventAction->AddEscapedSecondary(ekin);
      }
    }
  }
}
