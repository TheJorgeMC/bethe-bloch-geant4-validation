// ============================================================================
// Run.hh
// Subclase de G4Run que acumula, de forma segura entre hilos (cada hilo de
// trabajo posee su propia instancia local; G4 llama a Merge() sobre la del
// master), las sumas necesarias para las estadisticas de fin de run:
// media y desviacion estandar del deposito de energia (straggling), balance
// de energia por canales, y la comparacion contra G4EmCalculator al estilo
// del reporte de TestEm0.
// ============================================================================
#ifndef Run_hh
#define Run_hh 1

#include "G4Run.hh"
#include "globals.hh"

class DetectorConstruction;
class G4ParticleDefinition;

class Run : public G4Run
{
 public:
  explicit Run(const DetectorConstruction* det);
  ~Run() override = default;

  // Registrado por RunAction::BeginOfRunAction en los hilos de trabajo,
  // a partir del G4ParticleGun (patron de los ejemplos TestEm).
  void SetPrimary(const G4ParticleDefinition* particle, G4double energy);

  // Registrado por EventAction::EndOfEventAction, una vez por evento.
  //   exitEnergy < 0 significa que el primario NO salio del slab.
  void AddEvent(G4double edepPrimary, G4double edepSecondary,
                G4double escapedSecondary, G4double exitEnergy,
                G4double trackLength);

  void Merge(const G4Run* aRun) override;

  // Reporte de consola (solo lo llama el master en EndOfRunAction).
  void EndOfRun() const;

 private:
  const DetectorConstruction* fDetector;

  const G4ParticleDefinition* fParticle = nullptr;
  G4double fEkin = 0.;

  G4double fSumEdepTotal = 0.;   // primario + secundarios (por evento)
  G4double fSumEdepTotal2 = 0.;  // suma de cuadrados, para el straggling
  G4double fSumEdepPrimary = 0.;
  G4double fSumEdepSecondary = 0.;
  G4double fSumEscapedSecondary = 0.;
  G4double fSumExitEnergy = 0.;  // solo eventos en que el primario sale
  G4long fNExit = 0;             // numero de eventos con primario saliente
  G4double fSumTrackLength = 0.;
  G4long fNEvents = 0;
};

#endif
