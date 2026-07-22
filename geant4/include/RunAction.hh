// ============================================================================
// RunAction.hh
// Abre/cierra el archivo de G4AnalysisManager (CSV por defecto; ROOT si el
// nombre de archivo lleva extension .root), define el ntuple por evento y el
// histograma de dosis en profundidad, y delega el reporte de consola de fin
// de run en Run::EndOfRun() (solo en el master).
// ============================================================================
#ifndef RunAction_hh
#define RunAction_hh 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;
class Run;

class RunAction : public G4UserRunAction
{
 public:
  // primary == nullptr en el hilo master (patron TestEm: el master no tiene
  // generador de primarios en modo multi-hilo).
  RunAction(const DetectorConstruction* det,
            const PrimaryGeneratorAction* primary);
  ~RunAction() override = default;

  G4Run* GenerateRun() override;
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  void Book();  // definicion unica de ntuple e histogramas

  const DetectorConstruction* fDetector;
  const PrimaryGeneratorAction* fPrimary;
  Run* fRun = nullptr;
};

#endif
