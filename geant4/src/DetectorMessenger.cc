// ============================================================================
// DetectorMessenger.cc
// ============================================================================
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
  : G4UImessenger(), fDetector(det)
{
  fDirectory = new G4UIdirectory("/absorber/");
  fDirectory->SetGuidance("Absorber slab configuration");

  fMaterialCmd = new G4UIcmdWithAString("/absorber/material", this);
  fMaterialCmd->SetGuidance("NIST material of the slab (e.g. G4_WATER, G4_Al)");
  fMaterialCmd->SetParameterName("nistName", false);
  fMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fThicknessCmd = new G4UIcmdWithADoubleAndUnit("/absorber/thickness", this);
  fThicknessCmd->SetGuidance("Slab thickness along the beam axis (z)");
  fThicknessCmd->SetParameterName("thickness", false);
  fThicknessCmd->SetUnitCategory("Length");
  fThicknessCmd->SetRange("thickness>0.");
  fThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSizeXYCmd = new G4UIcmdWithADoubleAndUnit("/absorber/sizeXY", this);
  fSizeXYCmd->SetGuidance(
      "Side of the square transverse cross section of the slab. Must be "
      "large compared with the lateral spread to avoid energy leakage.");
  fSizeXYCmd->SetParameterName("sizeXY", false);
  fSizeXYCmd->SetUnitCategory("Length");
  fSizeXYCmd->SetRange("sizeXY>0.");
  fSizeXYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLayersCmd = new G4UIcmdWithAnInteger("/absorber/numberOfLayers", this);
  fLayersCmd->SetGuidance(
      "Number of layers (G4PVReplica along z). 1 = unsegmented slab; N>1 "
      "enables per-layer energy-deposit scoring (Bragg curve).");
  fLayersCmd->SetParameterName("N", false);
  fLayersCmd->SetRange("N>=1");
  fLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterialCmd;
  delete fThicknessCmd;
  delete fSizeXYCmd;
  delete fLayersCmd;
  delete fDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fMaterialCmd) {
    fDetector->SetAbsorberMaterial(newValue);
  } else if (command == fThicknessCmd) {
    fDetector->SetThickness(fThicknessCmd->GetNewDoubleValue(newValue));
  } else if (command == fSizeXYCmd) {
    fDetector->SetSizeXY(fSizeXYCmd->GetNewDoubleValue(newValue));
  } else if (command == fLayersCmd) {
    fDetector->SetNumberOfLayers(fLayersCmd->GetNewIntValue(newValue));
  } else {
    return;
  }

  // Flag the geometry as modified so that macro-driven changes take effect
  // without recompiling. In Geant4 11.x we use ReinitializeGeometry(), which
  // both marks the geometry as modified (equivalent to
  // GeometryHasBeenModified()) and forces Construct() to be re-run at the
  // next /run/beamOn or /run/initialize — required here because thickness/
  // layer changes rebuild solids, not just re-optimize navigation.
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
