// ============================================================================
// PrimaryGeneratorAction.hh
// G4ParticleGun with one monoenergetic proton per event, directed along +z,
// born just upstream of the slab entrance face. The energy is set per macro
// with the standard /gun/energy command (no additional messenger).
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
