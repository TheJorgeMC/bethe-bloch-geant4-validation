// ============================================================================
// EventAction.hh
// Acumuladores por evento, alimentados por SteppingAction. Al final del
// evento llena el ntuple, el histograma de dosis en profundidad y las sumas
// del objeto Run.
// ============================================================================
#ifndef EventAction_hh
#define EventAction_hh 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

class DetectorConstruction;

class EventAction : public G4UserEventAction
{
 public:
  explicit EventAction(const DetectorConstruction* det);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  // --- Interfaz para SteppingAction -------------------------------------
  void AddPrimaryEdep(G4double edep, G4int layer);
  void AddSecondaryEdep(G4double edep, G4int layer);
  void AddEscapedSecondary(G4double ekin) { fEscapedSecondary += ekin; }
  void AddPrimaryTrackLength(G4double len) { fPrimaryTrackLength += len; }
  void SetPrimaryExitEnergy(G4double ekin) { fPrimaryExitEnergy = ekin; }

 private:
  void AddLayerEdep(G4double edep, G4int layer);

  const DetectorConstruction* fDetector;

  G4double fEdepPrimary = 0.;
  G4double fEdepSecondary = 0.;
  G4double fEscapedSecondary = 0.;
  G4double fPrimaryExitEnergy = -1.;  // -1 = el primario no salio del slab
  G4double fPrimaryTrackLength = 0.;
  std::vector<G4double> fLayerEdep;   // deposito por capa (curva de Bragg)
};

#endif
