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
// El proton nace en z = -kGunOffsetFactor * espesor: entre la cara de
// entrada del slab (z = -0.5*t) y el borde del Envelope (z = -0.6*t),
// dentro del World de aire. El trayecto de aire recorrido antes de entrar
// al slab es de 0.15*t (p.ej. 0.75 mm para t = 5 mm), cuya perdida de
// energia es despreciable frente a las energias de este estudio.
constexpr G4double kGunOffsetFactor = 0.65;

// Energia por defecto (sobreescribible con /gun/energy en las macros).
constexpr G4double kDefaultEnergy = 150.0 * CLHEP::MeV;

// Un proton primario por evento (requisito de diseño): la estadistica se
// acumula con /run/beamOn N.
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
  // La posicion se recalcula en cada evento a partir del espesor actual del
  // slab, de modo que los cambios de geometria por macro (/absorber/thickness)
  // reposicionan automaticamente el vertice primario sin recompilar.
  const G4double z0 = -kGunOffsetFactor * fDetector->GetThickness();
  fGun->SetParticlePosition(G4ThreeVector(0., 0., z0));
  fGun->GeneratePrimaryVertex(event);
}
