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

#include "DetectorConstruction.hh"  // for kRegionName

namespace
{
// Default production cut: 1 mm (the historical Geant4 default). In water it
// corresponds to an electron threshold of ~350 keV: delta rays below that
// threshold are NOT generated explicitly and their energy is integrated
// into the continuous energy loss (restricted dE/dx) of the primary proton.
constexpr G4double kDefaultAbsorberCut = 1.0 * CLHEP::mm;
}  // namespace

PhysicsList::PhysicsList(G4int emOption)
  : G4VModularPhysicsList(), fAbsorberCut(kDefaultAbsorberCut)
{
  SetVerboseLevel(1);

  // --------------------------------------------------------------------------
  // Electromagnetic physics. option4 is the most accurate EM variant
  // available (low-energy Livermore/Penelope-type models and finer precision
  // parameters); option3 is offered as a model-sensitivity check. Any other
  // value falls back to option4 with a warning.
  // --------------------------------------------------------------------------
  if (emOption == 3) {
    RegisterPhysics(new G4EmStandardPhysics_option3());
    G4cout << "### PhysicsList: G4EmStandardPhysics_option3" << G4endl;
  } else {
    if (emOption != 4) {
      G4cout << "### PhysicsList: EM option '" << emOption
             << "' not recognized; using option4" << G4endl;
    }
    RegisterPhysics(new G4EmStandardPhysics_option4());
    G4cout << "### PhysicsList: G4EmStandardPhysics_option4" << G4endl;
  }

  // NOTE on G4EmExtraPhysics: deliberately OMITTED. That constructor adds
  // gamma-nuclear, muon-nuclear and synchrotron processes, none of which is
  // relevant to the electromagnetic energy loss of protons in a slab. The
  // Coulomb multiple scattering we DO need is already included in
  // G4EmStandardPhysics_option3/4 and does not require G4EmExtraPhysics.

  // Build the CSDA range table so that G4EmCalculator::GetCSDARange()
  // returns valid values in the end-of-run report (Run::EndOfRun).
  G4EmParameters::Instance()->SetBuildCSDARange(true);

  // Global default cut (World and everything not covered by the region).
  defaultCutValue = kDefaultAbsorberCut;

  // Messenger: /physics/absorberCut <value> <unit>
  fMessenger = new G4GenericMessenger(this, "/physics/",
                                      "Physics list parameters");
  fMessenger
      ->DeclareMethodWithUnit("absorberCut", "mm",
                              &PhysicsList::SetAbsorberCut,
                              "Production cut (length) for the AbsorberRegion "
                              "region. Equivalent to the native command "
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
  // If we are already initialized (Idle state) the kernel must be told to
  // rebuild the physics tables with the new cut at the next /run/beamOn;
  // in PreInit this is unnecessary (SetCuts() has not run yet).
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_Idle) {
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void PhysicsList::SetCuts()
{
  // ==========================================================================
  // PHYSICS PURPOSE (design, section 5): the production cut is expressed as
  // a LENGTH that Geant4 translates, material by material, into a threshold
  // energy T_cut through the range-energy relation. Collisions transferring
  // less than T_cut do not spawn an explicit delta ray: their energy is
  // integrated continuously into the proton step (RESTRICTED dE/dx).
  // Collisions transferring more than T_cut spawn a secondary electron that
  // is transported and deposits its energy elsewhere. Therefore:
  //   dE/dx(unrestricted) = dE/dx(restricted, T_cut)
  //                        + energy carried by explicit delta rays
  // and the value reported by the continuous process DEPENDS on the chosen
  // cut: it is a simulation design choice, not a physical constant.
  // ==========================================================================
  SetCutsWithDefault();

  G4Region* region = G4RegionStore::GetInstance()->GetRegion(
      DetectorConstruction::kRegionName, false);
  if (region != nullptr) {
    auto* cuts = new G4ProductionCuts();
    // Same length cut for all relevant secondary species. For the delta-ray
    // study the decisive one is "e-".
    cuts->SetProductionCut(fAbsorberCut, "gamma");
    cuts->SetProductionCut(fAbsorberCut, "e-");
    cuts->SetProductionCut(fAbsorberCut, "e+");
    cuts->SetProductionCut(fAbsorberCut, "proton");
    region->SetProductionCuts(cuts);

    G4cout << "### PhysicsList: production cut in "
           << DetectorConstruction::kRegionName << " = "
           << G4BestUnit(fAbsorberCut, "Length") << G4endl;
  } else {
    G4cout << "### PhysicsList: warning — region "
           << DetectorConstruction::kRegionName
           << " does not exist yet; only the default cut was applied"
           << G4endl;
  }

  if (verboseLevel > 0) DumpCutValuesTable();
}
