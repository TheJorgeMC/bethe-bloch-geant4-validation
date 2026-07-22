// ============================================================================
// PrimaryGeneratorAction.cc
// ============================================================================
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

namespace
{
// The proton is born at z = -kGunOffsetFactor * thickness: between the slab
// entrance face (z = -0.5*t) and the Envelope edge (z = -0.6*t), inside the
// air World. The air path traversed before entering the slab is 0.15*t
// (e.g. 0.75 mm for t = 5 mm), whose energy loss is negligible compared
// with the energies of this study.
constexpr G4double kGunOffsetFactor = 0.65;

// Default energy (overridable with /gun/energy in the macros).
constexpr G4double kDefaultEnergy = 150.0 * CLHEP::MeV;

// One primary proton per event (design requirement): statistics are
// accumulated with /run/beamOn N.
constexpr G4int kParticlesPerEvent = 1;
}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction(const DetectorConstruction* det)
  : G4VUserPrimaryGeneratorAction(), fDetector(det)
{
  fGun = new G4ParticleGun(kParticlesPerEvent);

  G4ParticleDefinition* proton =
      G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGun->SetParticleDefinition(proton);
  fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fGun->SetParticleEnergy(kDefaultEnergy);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // The position is recomputed every event from the current slab thickness,
  // so macro-driven geometry changes (/absorber/thickness) automatically
  // reposition the primary vertex without recompiling.
  const G4double z0 = -kGunOffsetFactor * fDetector->GetThickness();
  fGun->SetParticlePosition(G4ThreeVector(0., 0., z0));
  fGun->GeneratePrimaryVertex(event);
}
