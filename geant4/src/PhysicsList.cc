// ============================================================================
// PhysicsList.cc
// ============================================================================
#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmParameters.hh"
#include "G4GenericMessenger.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"  // para kRegionName

namespace
{
// Corte de produccion por defecto: 1 mm (el default historico de Geant4).
// En agua equivale a un umbral de ~350 keV para electrones: rayos delta por
// debajo de ese umbral NO se generan explicitamente y su energia se integra
// en la perdida continua (dE/dx restringido) del proton primario.
constexpr G4double kDefaultAbsorberCut = 1.0 * CLHEP::mm;
}  // namespace

PhysicsList::PhysicsList(G4int emOption)
  : G4VModularPhysicsList(), fAbsorberCut(kDefaultAbsorberCut)
{
  SetVerboseLevel(1);

  // --------------------------------------------------------------------------
  // Fisica electromagnetica. option4 es la variante EM mas precisa disponible
  // (modelos de baja energia tipo Livermore/Penelope y parametros de precision
  // mas finos); option3 se ofrece como chequeo de sensibilidad al modelo.
  // Cualquier otro valor cae en option4 con un aviso.
  // --------------------------------------------------------------------------
  if (emOption == 3) {
    RegisterPhysics(new G4EmStandardPhysics_option3());
    G4cout << "### PhysicsList: G4EmStandardPhysics_option3" << G4endl;
  } else {
    if (emOption != 4) {
      G4cout << "### PhysicsList: opcion EM '" << emOption
             << "' no reconocida; se usa option4" << G4endl;
    }
    RegisterPhysics(new G4EmStandardPhysics_option4());
    G4cout << "### PhysicsList: G4EmStandardPhysics_option4" << G4endl;
  }

  // NOTA sobre G4EmExtraPhysics: se OMITE deliberadamente. Ese constructor
  // añade procesos gamma-nucleares, muon-nucleares y de synchrotron, ninguno
  // relevante para la perdida de energia electromagnetica de protones en un
  // slab. La dispersion multiple de Coulomb (que SI necesitamos) ya viene
  // incluida en G4EmStandardPhysics_option3/4; no requiere G4EmExtraPhysics.

  // Construir la tabla de rango CSDA para que G4EmCalculator::GetCSDARange()
  // devuelva valores validos en el reporte de fin de run (Run::EndOfRun).
  G4EmParameters::Instance()->SetBuildCSDARange(true);

  // Corte por defecto global (World y todo lo no cubierto por la region).
  defaultCutValue = kDefaultAbsorberCut;

  // Messenger: /physics/absorberCut <valor> <unidad>
  fMessenger = new G4GenericMessenger(this, "/physics/",
                                      "Parametros de la physics list");
  fMessenger
      ->DeclareMethodWithUnit("absorberCut", "mm",
                              &PhysicsList::SetAbsorberCut,
                              "Corte de produccion (longitud) para la region "
                              "AbsorberRegion. Equivalente al comando nativo "
                              "/run/setCutForRegion AbsorberRegion <v> <u>.")
      .SetStates(G4State_PreInit, G4State_Idle);
}

PhysicsList::~PhysicsList()
{
  delete fMessenger;
}

void PhysicsList::SetAbsorberCut(G4double value)
{
  fAbsorberCut = value;
  // Si ya estamos inicializados (estado Idle), hay que avisar al kernel para
  // que reconstruya las tablas fisicas con el nuevo corte en el siguiente
  // /run/beamOn; en PreInit no hace falta (SetCuts() aun no ha corrido).
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Idle) {
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void PhysicsList::SetCuts()
{
  // ==========================================================================
  // PROPOSITO FISICO (diseño, seccion 5): el corte de produccion se expresa
  // como una LONGITUD que Geant4 traduce, material por material, a una
  // energia umbral T_cut via la relacion rango-energia. Colisiones con
  // transferencia < T_cut no generan rayo delta explicito: su energia se
  // integra de forma continua en el paso del proton (dE/dx RESTRINGIDO).
  // Colisiones con transferencia > T_cut generan un electron secundario que
  // se transporta y deposita su energia en otro lugar. Por tanto:
  //   dE/dx(no restringido) = dE/dx(restringido, T_cut)
  //                          + energia llevada por deltas explicitos
  // y el valor reportado por el proceso continuo DEPENDE del corte elegido:
  // es una decision de diseño de la simulacion, no una constante fisica.
  // ==========================================================================
  SetCutsWithDefault();

  G4Region* region = G4RegionStore::GetInstance()->GetRegion(
      DetectorConstruction::kRegionName, false);
  if (region != nullptr) {
    auto* cuts = new G4ProductionCuts();
    // El mismo corte de longitud para todas las especies secundarias
    // relevantes. Para el estudio de rayos delta el decisivo es "e-".
    cuts->SetProductionCut(fAbsorberCut, "gamma");
    cuts->SetProductionCut(fAbsorberCut, "e-");
    cuts->SetProductionCut(fAbsorberCut, "e+");
    cuts->SetProductionCut(fAbsorberCut, "proton");
    region->SetProductionCuts(cuts);

    G4cout << "### PhysicsList: corte de produccion en "
           << DetectorConstruction::kRegionName << " = "
           << G4BestUnit(fAbsorberCut, "Length") << G4endl;
  } else {
    G4cout << "### PhysicsList: aviso — la region "
           << DetectorConstruction::kRegionName
           << " no existe todavia; solo se aplico el corte por defecto"
           << G4endl;
  }

  if (verboseLevel > 0) DumpCutValuesTable();
}
