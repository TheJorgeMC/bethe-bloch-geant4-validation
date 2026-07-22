// ============================================================================
// EventAction.cc
// ============================================================================
#include "EventAction.hh"
#include "AnalysisDefs.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(const DetectorConstruction* det)
  : G4UserEventAction(), fDetector(det)
{
}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdepPrimary = 0.;
  fEdepSecondary = 0.;
  fEscapedSecondary = 0.;
  fPrimaryExitEnergy = -1.;
  fPrimaryTrackLength = 0.;
  fLayerEdep.assign(fDetector->GetNumberOfLayers(), 0.);
}

void EventAction::AddLayerEdep(G4double edep, G4int layer)
{
  if (layer >= 0 && layer < static_cast<G4int>(fLayerEdep.size())) {
    fLayerEdep[layer] += edep;
  }
}

void EventAction::AddPrimaryEdep(G4double edep, G4int layer)
{
  fEdepPrimary += edep;
  AddLayerEdep(edep, layer);
}

void EventAction::AddSecondaryEdep(G4double edep, G4int layer)
{
  fEdepSecondary += edep;
  AddLayerEdep(edep, layer);
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Energia cinetica del primario al nacer, leida del vertice del evento
  // (robusto frente a cambios de /gun/energy entre runs).
  G4double eIncident = 0.;
  if (event->GetPrimaryVertex() != nullptr &&
      event->GetPrimaryVertex()->GetPrimary() != nullptr) {
    eIncident = event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
  }

  // --- Ntuple: una fila por evento, en MeV y mm -------------------------
  auto* am = G4AnalysisManager::Instance();
  am->FillNtupleDColumn(analysis::kColEIncident, eIncident / MeV);
  am->FillNtupleDColumn(analysis::kColEdepPrimary, fEdepPrimary / MeV);
  am->FillNtupleDColumn(analysis::kColEdepSecondary, fEdepSecondary / MeV);
  am->FillNtupleDColumn(analysis::kColEEscapedSec, fEscapedSecondary / MeV);
  am->FillNtupleDColumn(analysis::kColEExitPrimary,
                        (fPrimaryExitEnergy >= 0.) ? fPrimaryExitEnergy / MeV
                                                   : -1.);
  am->FillNtupleDColumn(analysis::kColTrackLenPrimary,
                        fPrimaryTrackLength / mm);
  am->AddNtupleRow();

  // --- Histograma de dosis en profundidad (suma sobre eventos) ----------
  // Cada capa aporta su deposito en el centro geometrico del bin. Con
  // numberOfLayers = 1 degenera en un unico bin (deposito total).
  const G4int nLayers = static_cast<G4int>(fLayerEdep.size());
  const G4double layerThicknessMM =
      (fDetector->GetThickness() / mm) / nLayers;
  for (G4int i = 0; i < nLayers; ++i) {
    if (fLayerEdep[i] > 0.) {
      const G4double zCenter = (i + 0.5) * layerThicknessMM;
      am->FillH1(analysis::kH1Depth, zCenter, fLayerEdep[i] / MeV);
    }
  }

  // --- Sumas del run (fusionadas entre hilos en Run::Merge) --------------
  auto* run = static_cast<Run*>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->AddEvent(fEdepPrimary, fEdepSecondary, fEscapedSecondary,
                fPrimaryExitEnergy, fPrimaryTrackLength);
}
