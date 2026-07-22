// ============================================================================
// DetectorConstruction.cc
// ============================================================================
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"

namespace
{
// Size factors relative to the slab (no scattered "magic numbers").
// The World is 1.4x the slab along each axis; the Envelope, 1.2x. The
// Envelope is an optional hook for scoring overlays or for changing the
// background material without touching the World (design document, sec. 2.1).
constexpr G4double kWorldScale = 1.4;
constexpr G4double kEnvelopeScale = 1.2;

// Geometry defaults
constexpr G4double kDefaultThickness = 5.0 * CLHEP::mm;
constexpr G4double kDefaultSizeXY = 10.0 * CLHEP::cm;
constexpr G4int kDefaultNumberOfLayers = 1;
const G4String kDefaultMaterial = "G4_WATER";

// Overlap check in every G4PVPlacement (design requirement).
constexpr G4bool kCheckOverlaps = true;
}  // namespace

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fThickness(kDefaultThickness),
    fSizeXY(kDefaultSizeXY),
    fNumberOfLayers(kDefaultNumberOfLayers)
{
  fAbsorberMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial(kDefaultMaterial);
  fMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

void DetectorConstruction::CleanupGeometry()
{
  // When rebuilding the geometry after a macro-driven change the stores must
  // be cleaned; otherwise stale volumes are left dangling and Geant4 aborts.
  // Before cleaning, unregister the root logical volume from the region so
  // that "AbsorberRegion" is not left holding a dangling pointer.
  if (fRegion != nullptr && fAbsorberLV != nullptr) {
    fRegion->RemoveRootLogicalVolume(fAbsorberLV, false);
  }
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  fWorldLV = nullptr;
  fWorldPV = nullptr;
  fEnvelopeLV = nullptr;
  fAbsorberLV = nullptr;
  fLayerLV = nullptr;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (fWorldPV != nullptr) CleanupGeometry();

  auto* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

  // --------------------------------------------------------------------------
  // World: air box, 1.4x the slab. Physics note: the primary proton is born
  // in air at ~0.15*thickness upstream of the entrance face (see
  // PrimaryGeneratorAction); the energy loss along that air path is
  // negligible (<1 keV for the energies of this study). To remove it
  // entirely, change the World material to G4_Galactic here.
  // --------------------------------------------------------------------------
  const G4double worldHalfXY = 0.5 * kWorldScale * fSizeXY;
  const G4double worldHalfZ = 0.5 * kWorldScale * fThickness;
  auto* worldSolid = new G4Box("World", worldHalfXY, worldHalfXY, worldHalfZ);
  fWorldLV = new G4LogicalVolume(worldSolid, air, "World");
  fWorldPV = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLV, "World", nullptr,
                               false, 0, kCheckOverlaps);

  // --------------------------------------------------------------------------
  // Envelope: same material as the World, 1.2x the slab.
  // --------------------------------------------------------------------------
  const G4double envHalfXY = 0.5 * kEnvelopeScale * fSizeXY;
  const G4double envHalfZ = 0.5 * kEnvelopeScale * fThickness;
  auto* envSolid = new G4Box("Envelope", envHalfXY, envHalfXY, envHalfZ);
  fEnvelopeLV = new G4LogicalVolume(envSolid, air, "Envelope");
  new G4PVPlacement(nullptr, G4ThreeVector(), fEnvelopeLV, "Envelope", fWorldLV,
                    false, 0, kCheckOverlaps);

  // --------------------------------------------------------------------------
  // Absorber slab: G4Box with the beam axis along z.
  // If fNumberOfLayers > 1, the slab container is subdivided into equal
  // layers with G4PVReplica along kZAxis; each layer carries a copy number
  // (copyNo) that SteppingAction uses for depth scoring.
  // --------------------------------------------------------------------------
  auto* slabSolid =
      new G4Box("Absorber", 0.5 * fSizeXY, 0.5 * fSizeXY, 0.5 * fThickness);
  fAbsorberLV = new G4LogicalVolume(slabSolid, fAbsorberMaterial, "Absorber");
  new G4PVPlacement(nullptr, G4ThreeVector(), fAbsorberLV, "Absorber",
                    fEnvelopeLV, false, 0, kCheckOverlaps);

  if (fNumberOfLayers > 1) {
    const G4double layerThickness = fThickness / fNumberOfLayers;
    auto* layerSolid = new G4Box("Layer", 0.5 * fSizeXY, 0.5 * fSizeXY,
                                 0.5 * layerThickness);
    // The layer must use the same material as the container (a replica
    // fills its mother volume completely).
    fLayerLV = new G4LogicalVolume(layerSolid, fAbsorberMaterial, "Layer");
    // G4PVReplica takes no overlap-check flag: the regular subdivision is
    // exact by construction.
    new G4PVReplica("Layer", fLayerLV, fAbsorberLV, kZAxis, fNumberOfLayers,
                    layerThickness);
  } else {
    fLayerLV = nullptr;
  }

  // --------------------------------------------------------------------------
  // Dedicated region for the slab, enabling local production cuts
  // (/physics/absorberCut, see PhysicsList::SetCuts()).
  // --------------------------------------------------------------------------
  fRegion = G4RegionStore::GetInstance()->GetRegion(kRegionName, false);
  if (fRegion == nullptr) {
    fRegion = new G4Region(kRegionName);
  }
  fRegion->AddRootLogicalVolume(fAbsorberLV);

  // --- Minimal visualization attributes ---
  fWorldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* envVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.1));
  envVis->SetForceWireframe(true);
  fEnvelopeLV->SetVisAttributes(envVis);
  auto* slabVis = new G4VisAttributes(G4Colour(0.2, 0.4, 0.8, 0.3));
  slabVis->SetForceSolid(true);
  fAbsorberLV->SetVisAttributes(slabVis);

  G4cout << "### DetectorConstruction: slab of "
         << fAbsorberMaterial->GetName() << ", thickness "
         << G4BestUnit(fThickness, "Length") << ", cross section "
         << G4BestUnit(fSizeXY, "Length") << " x "
         << G4BestUnit(fSizeXY, "Length") << ", " << fNumberOfLayers
         << " layer(s)" << G4endl;

  return fWorldPV;
}

// ----------------------------------------------------------------------------
// Messenger setters. The actual rebuild is triggered by DetectorMessenger
// after calling these methods (see DetectorMessenger.cc).
// ----------------------------------------------------------------------------
void DetectorConstruction::SetAbsorberMaterial(const G4String& nistName)
{
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(nistName);
  if (mat == nullptr) {
    G4cerr << "### DetectorConstruction: unknown NIST material '"
           << nistName << "'. Keeping "
           << fAbsorberMaterial->GetName() << G4endl;
    return;
  }
  fAbsorberMaterial = mat;
}

void DetectorConstruction::SetThickness(G4double value)
{
  if (value <= 0.) {
    G4cerr << "### DetectorConstruction: invalid thickness, ignored" << G4endl;
    return;
  }
  fThickness = value;
}

void DetectorConstruction::SetSizeXY(G4double value)
{
  if (value <= 0.) {
    G4cerr << "### DetectorConstruction: invalid sizeXY, ignored" << G4endl;
    return;
  }
  fSizeXY = value;
}

void DetectorConstruction::SetNumberOfLayers(G4int n)
{
  if (n < 1) {
    G4cerr << "### DetectorConstruction: numberOfLayers < 1, ignored"
           << G4endl;
    return;
  }
  fNumberOfLayers = n;
}
