// ============================================================================
// slab.cc — main()
// Validacion de la ecuacion de Bethe-Bloch: protones monoenergeticos sobre un
// slab absorbedor homogeneo (diseño Fase 5).
//
// Uso:
//   ./slab                      -> sesion interactiva (ejecuta macros/vis.mac)
//   ./slab macros/vis.mac       -> sesion interactiva con esa macro (los
//                                  drivers Qt/OpenGL requieren una sesion
//                                  G4UIExecutive creada ANTES de /vis/open;
//                                  cualquier macro cuyo nombre contenga "vis"
//                                  se ejecuta por eso en modo interactivo)
//   ./slab macros/run.mac       -> modo batch, EM option4 (default)
//   ./slab macros/run.mac 3     -> modo batch, EM option3 (chequeo de
//                                  sensibilidad al modelo EM)
//
// G4RunManagerFactory crea automaticamente un run manager multi-hilo si
// Geant4 fue compilado con soporte MT; en caso contrario cae a secuencial.
// El numero de hilos puede fijarse con /run/numberOfThreads en la macro
// (antes de /run/initialize) o con la variable G4FORCENUMBEROFTHREADS.
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

int main(int argc, char** argv)
{
  // --- Argumentos: [macro] [emOption] ---
  G4String macroFile;
  G4int emOption = 4;  // default: G4EmStandardPhysics_option4
  if (argc > 1) macroFile = argv[1];
  if (argc > 2) emOption = std::atoi(argv[2]);

  // Sesion interactiva si no se paso macro O si la macro es de visualizacion
  // (nombre conteniendo "vis"): los drivers interactivos Qt/OpenGL exigen que
  // la sesion exista antes de crear el sistema grafico con /vis/open.
  const G4bool isVisMacro =
      !macroFile.empty() && macroFile.find("vis") != std::string::npos;
  G4UIExecutive* ui = nullptr;
  if (macroFile.empty() || isVisMacro) {
    ui = new G4UIExecutive(argc, argv);
  }

  // MT si esta disponible; secuencial como fallback (G4RunManagerType::Default).
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  auto* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList(emOption));
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // Visualizacion (necesaria para vis.mac; inocua en batch).
  auto* visManager = new G4VisExecutive("quiet");
  visManager->Initialize();

  G4UImanager* uiManager = G4UImanager::GetUIpointer();

  if (ui == nullptr) {
    // Modo batch: ejecutar la macro indicada.
    uiManager->ApplyCommand("/control/execute " + macroFile);
  } else {
    // Modo interactivo: ejecutar la macro indicada (o macros/vis.mac por
    // defecto) con la sesion ya creada, y abrir la sesion.
    const G4String visMacro = isVisMacro ? macroFile : G4String("macros/vis.mac");
    uiManager->ApplyCommand("/control/execute " + visMacro);
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
  return 0;
}
