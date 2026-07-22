// ============================================================================
// SteppingAction.hh
// Scoring manual paso a paso dentro del slab: separa la energia depositada
// por el proton primario (ParentID == 0) de la de los secundarios (rayos
// delta y su progenie), contabiliza la energia que los secundarios sacan del
// slab, y registra energia de salida y longitud de traza del primario.
// Este es el mecanismo que permite reconstruir el balance
//   dE/dx(no restringido) = dE/dx(restringido) + energia de deltas explicitos
// del documento de diseño (seccion 5.3).
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
