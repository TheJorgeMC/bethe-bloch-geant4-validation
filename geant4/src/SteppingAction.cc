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
  // Solo nos interesan los pasos que OCURREN dentro del slab: el volumen del
  // punto pre-step es el que define donde se deposito la energia del paso.
  const G4StepPoint* prePoint = step->GetPreStepPoint();
  const G4LogicalVolume* preLV =
      prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (!fDetector->IsAbsorber(preLV)) return;

  G4Track* track = step->GetTrack();
  const G4bool isPrimary = (track->GetParentID() == 0);

  // Indice de capa para el scoring por profundidad. Con G4PVReplica el
  // numero de replica es el copy number del touchable en profundidad 0.
  // Con slab sin segmentar (numberOfLayers == 1) el copy number es 0.
  const G4int layer = (fDetector->GetNumberOfLayers() > 1)
                          ? prePoint->GetTouchableHandle()->GetCopyNumber()
                          : 0;

  // --------------------------------------------------------------------------
  // 1) Energia depositada en este paso.
  // FISICA: para el primario, GetTotalEnergyDeposit() contiene la perdida
  // CONTINUA del proceso de ionizacion — es decir, el dE/dx RESTRINGIDO con
  // T_cut fijado por el corte de produccion de la region — mas cualquier
  // deposito discreto bajo umbral. La energia cedida a rayos delta EXPLICITOS
  // NO esta incluida aqui: viaja con la traza secundaria y aparecera en el
  // contador de secundarios (si el delta se frena dentro del slab) o en el
  // de energia escapada (si sale).
  // --------------------------------------------------------------------------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0.) {
    if (isPrimary) {
      fEventAction->AddPrimaryEdep(edep, layer);
    } else {
      fEventAction->AddSecondaryEdep(edep, layer);
    }
  }

  // 2) Longitud de traza del primario dentro del slab (para normalizar por
  // camino real en lugar de espesor geometrico si la dispersion multiple es
  // apreciable).
  if (isPrimary) {
    fEventAction->AddPrimaryTrackLength(step->GetStepLength());
  }

  // --------------------------------------------------------------------------
  // 3) Cruce de frontera hacia fuera del slab.
  // Si el post-step esta en una frontera geometrica y el siguiente volumen ya
  // no pertenece al slab (ojo: pasar de una capa replicada a la siguiente
  // sigue siendo slab, IsAbsorber lo cubre), la traza esta saliendo.
  // --------------------------------------------------------------------------
  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    const G4VPhysicalVolume* nextPV = track->GetNextVolume();
    const G4LogicalVolume* nextLV =
        (nextPV != nullptr) ? nextPV->GetLogicalVolume() : nullptr;
    if (nextLV == nullptr || !fDetector->IsAbsorber(nextLV)) {
      const G4double ekin = track->GetKineticEnergy();
      if (isPrimary) {
        // Energia cinetica del primario al salir del slab (por la cara
        // opuesta o, tras mucha dispersion, por un lateral: se registra
        // igual; el ntuple permite filtrar por longitud de traza).
        fEventAction->SetPrimaryExitEnergy(ekin);
      } else if (ekin > 0.) {
        // Energia que un secundario (tipicamente un rayo delta o un foton
        // de bremsstrahlung/fluorescencia) SACA del slab sin depositarla.
        // Esta es la componente que hay que sumar al deposito local para
        // recuperar el dE/dx no restringido (balance de la seccion 5.3).
        // NOTA: si la misma traza reentrara al slab y volviera a salir se
        // contaria dos veces; con esta geometria (slab rodeado de aire) el
        // efecto es despreciable y se acepta como aproximacion.
        fEventAction->AddEscapedSecondary(ekin);
      }
    }
  }
}
