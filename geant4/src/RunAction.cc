// ============================================================================
// RunAction.cc
// ============================================================================
#include "RunAction.hh"
#include "AnalysisDefs.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"

namespace
{
// Nombre de archivo por defecto; se sobreescribe con /analysis/setFileName.
const G4String kDefaultFileName = "slab_out";
// Tipo de archivo por defecto. CSV es directamente legible con pandas;
// para un unico archivo fusionado en modo MT, usar extension .root en
// /analysis/setFileName (p.ej. dedx_150MeV.root) — ver README.
const G4String kDefaultFileType = "csv";
}  // namespace

RunAction::RunAction(const DetectorConstruction* det,
                     const PrimaryGeneratorAction* primary)
  : G4UserRunAction(), fDetector(det), fPrimary(primary)
{
  Book();
}

void RunAction::Book()
{
  auto* am = G4AnalysisManager::Instance();
  am->SetDefaultFileType(kDefaultFileType);
  am->SetFileName(kDefaultFileName);
  am->SetVerboseLevel(0);
  // La fusion de ntuples entre hilos solo esta soportada para ROOT; para CSV
  // se generan archivos por hilo (_t0, _t1, ...) que se concatenan en el
  // post-analisis (ver README).
  am->SetNtupleMerging(true);

  // --- H1 kH1Depth: dosis en profundidad. Se crea con un binning provisional
  // y se re-binnea en BeginOfRunAction con el numero de capas y espesor
  // vigentes (SetH1 permite redefinir bins entre runs).
  const G4int h1 = am->CreateH1(
      "DepthEdep", "Energia depositada vs profundidad;z (mm);Edep (MeV)", 1,
      0., 1.);
  if (h1 != analysis::kH1Depth) {
    G4Exception("RunAction::Book()", "SlabAnalysis001", FatalException,
                "El id del histograma DepthEdep no coincide con kH1Depth");
  }

  // --- Ntuple por evento (una fila por proton primario) ---
  am->CreateNtuple("slab", "Balance de energia por evento (MeV, mm)");
  am->CreateNtupleDColumn("E_incidente");        // analysis::kColEIncident
  am->CreateNtupleDColumn("E_deposit_primario"); // analysis::kColEdepPrimary
  am->CreateNtupleDColumn("E_deposit_secundarios");
  am->CreateNtupleDColumn("E_transportada_fuera");
  am->CreateNtupleDColumn("E_salida_primario");
  am->CreateNtupleDColumn("longitud_traza");
  am->FinishNtuple();
}

G4Run* RunAction::GenerateRun()
{
  fRun = new Run(fDetector);
  return fRun;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  // En los hilos de trabajo, registrar particula y energia del gun en el Run
  // local; Merge() las propaga al Run del master (patron TestEm0/TestEm1).
  if (fPrimary != nullptr && fRun != nullptr) {
    fRun->SetPrimary(fPrimary->GetGun()->GetParticleDefinition(),
                     fPrimary->GetGun()->GetParticleEnergy());
  }

  auto* am = G4AnalysisManager::Instance();

  // Re-binnear el histograma de profundidad con la geometria vigente:
  // un bin por capa si el slab esta segmentado; un unico bin si no.
  const G4int nLayers = fDetector->GetNumberOfLayers();
  const G4double thicknessMM = fDetector->GetThickness() / mm;
  am->SetH1(analysis::kH1Depth, nLayers, 0., thicknessMM);

  am->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* am = G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();

  // El reporte de consola con estadisticas fusionadas de todos los hilos
  // solo tiene sentido en el master.
  if (IsMaster() && fRun != nullptr) {
    fRun->EndOfRun();
  }
}
