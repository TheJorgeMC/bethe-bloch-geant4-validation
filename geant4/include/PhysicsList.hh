// ============================================================================
// PhysicsList.hh
// Lista modular minima: SOLO fisica electromagnetica estandar
// (G4EmStandardPhysics_option4 por defecto; option3 seleccionable como
// chequeo de sensibilidad via argumento de linea de comandos en slab.cc:
//   ./slab run.mac 3
// ).
//
// Deliberadamente SIN fisica hadronica: para protones atravesando un slab
// delgado, las reacciones nucleares inelasticas son un efecto de segundo
// orden para la validacion de dE/dx y solo añadirian un canal de perdida de
// primarios y tiempo de computo (decision documentada en el diseño, seccion 4).
//
// Corte de produccion especifico del slab:
//   /physics/absorberCut <valor> <unidad>   (default 1 mm)
// aplicado a la G4Region "AbsorberRegion" en SetCuts(). Este corte fija el
// umbral T_cut de produccion explicita de rayos delta y, por tanto, define
// que dE/dx restringido usa el proceso continuo de ionizacion (diseño,
// seccion 5).
// ============================================================================
#ifndef PhysicsList_hh
#define PhysicsList_hh 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4GenericMessenger;

class PhysicsList : public G4VModularPhysicsList
{
 public:
  // emOption: 3 -> G4EmStandardPhysics_option3, 4 -> option4 (default)
  explicit PhysicsList(G4int emOption = 4);
  ~PhysicsList() override;

  void SetCuts() override;

  // Invocado por el messenger /physics/absorberCut
  void SetAbsorberCut(G4double value);

 private:
  G4double fAbsorberCut;  // corte de produccion (longitud) para AbsorberRegion
  G4GenericMessenger* fMessenger = nullptr;
};

#endif
