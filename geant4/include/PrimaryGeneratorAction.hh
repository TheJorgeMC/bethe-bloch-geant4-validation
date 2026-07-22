// ============================================================================
// PrimaryGeneratorAction.hh
// G4ParticleGun con un proton monoenergetico por evento, dirigido en +z,
// nacido justo antes de la cara de entrada del slab. La energia se fija por
// macro con el comando estandar /gun/energy (sin messenger adicional).
// ============================================================================
#ifndef PrimaryGeneratorAction_hh
#define PrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  explicit PrimaryGeneratorAction(const DetectorConstruction* det);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  const G4ParticleGun* GetGun() const { return fGun; }

 private:
  G4ParticleGun* fGun = nullptr;
  const DetectorConstruction* fDetector;
};

#endif
