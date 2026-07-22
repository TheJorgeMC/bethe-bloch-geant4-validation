// ============================================================================
// slab.cc — main()
// Bethe-Bloch equation validation: monoenergetic protons on a homogeneous
// absorber slab (Phase 5 design).
//
// Usage:
//   ./slab                      -> interactive session (runs macros/vis.mac)
//   ./slab macros/vis.mac       -> interactive session with that macro (the
//                                  Qt/OpenGL drivers require a G4UIExecutive
//                                  session created BEFORE /vis/open; any
//                                  macro whose name contains "vis" is
//                                  therefore run in interactive mode)
//   ./slab macros/run.mac       -> batch mode, EM option4 (default)
//   ./slab macros/run.mac 3     -> batch mode, EM option3 (EM-model
//                                  sensitivity check)
//
// G4RunManagerFactory automatically creates a multi-threaded run manager if
// Geant4 was built with MT support; otherwise it falls back to sequential.
// The number of threads can be set with /run/numberOfThreads in the macro
// (before /run/initialize) or with the G4FORCENUMBEROFTHREADS environment
// variable.
// ============================================================================
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "globals.hh"

#include <cstdlib>
#include <string>

int main(int argc, char** argv)
{
  // --- Arguments: [macro] [emOption] ---
  G4String macroFile;
  G4int emOption = 4;  // default: G4EmStandardPhysics_option4
  if (argc > 1) macroFile = argv[1];
  if (argc > 2) emOption = std::atoi(argv[2]);

  // Interactive session if no macro was given OR if the macro is a
  // visualization one (name containing "vis"): the interactive Qt/OpenGL
  // drivers require the session to exist before the graphics system is
  // created with /vis/open.
  const G4bool isVisMacro =
      !macroFile.empty() && macroFile.find("vis") != std::string::npos;
  G4UIExecutive* ui = nullptr;
  if (macroFile.empty() || isVisMacro) {
    ui = new G4UIExecutive(argc, argv);
  }

  // MT if available; sequential as fallback (G4RunManagerType::Default).
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  auto* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList(emOption));
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // Visualization (needed for vis.mac; harmless in batch).
  auto* visManager = new G4VisExecutive("quiet");
  visManager->Initialize();

  G4UImanager* uiManager = G4UImanager::GetUIpointer();

  if (ui == nullptr) {
    // Batch mode: execute the given macro.
    uiManager->ApplyCommand("/control/execute " + macroFile);
  } else {
    // Interactive mode: execute the given macro (or macros/vis.mac by
    // default) with the session already created, then start the session.
    const G4String visMacro = isVisMacro ? macroFile : G4String("macros/vis.mac");
    uiManager->ApplyCommand("/control/execute " + visMacro);
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}
