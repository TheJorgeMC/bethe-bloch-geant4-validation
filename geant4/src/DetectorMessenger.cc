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
  fDirectory->SetGuidance("Configuracion del slab absorbedor");

  fMaterialCmd = new G4UIcmdWithAString("/absorber/material", this);
  fMaterialCmd->SetGuidance("Material NIST del slab (p.ej. G4_WATER, G4_Al)");
  fMaterialCmd->SetParameterName("nistName", false);
  fMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fThicknessCmd = new G4UIcmdWithADoubleAndUnit("/absorber/thickness", this);
  fThicknessCmd->SetGuidance("Espesor del slab a lo largo del eje del haz (z)");
  fThicknessCmd->SetParameterName("thickness", false);
  fThicknessCmd->SetUnitCategory("Length");
  fThicknessCmd->SetRange("thickness>0.");
  fThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSizeXYCmd = new G4UIcmdWithADoubleAndUnit("/absorber/sizeXY", this);
  fSizeXYCmd->SetGuidance(
      "Lado de la seccion transversal cuadrada del slab. Debe ser grande "
      "frente a la dispersion lateral para evitar fugas de energia.");
  fSizeXYCmd->SetParameterName("sizeXY", false);
  fSizeXYCmd->SetUnitCategory("Length");
  fSizeXYCmd->SetRange("sizeXY>0.");
  fSizeXYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLayersCmd = new G4UIcmdWithAnInteger("/absorber/numberOfLayers", this);
  fLayersCmd->SetGuidance(
      "Numero de capas (G4PVReplica en z). 1 = slab sin segmentar; N>1 "
      "activa el scoring de energia depositada por capa (curva de Bragg).");
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

  // Marca la geometría como modificada para que el cambio por macro surta
  // efecto sin recompilar. En Geant4 11.x usamos ReinitializeGeometry(),
  // que ademas de marcar la geometria como modificada (equivalente a
  // GeometryHasBeenModified()) fuerza a re-ejecutar Construct() en el
  // siguiente /run/beamOn o /run/initialize — necesario aqui porque los
  // cambios de espesor/capas requieren reconstruir los solidos, no solo
  // reoptimizar la navegacion.
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
