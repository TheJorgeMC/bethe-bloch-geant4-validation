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
// Factores de tamaño relativos al slab (sin "magic numbers" dispersos).
// El World es 1.4x el slab en cada eje; el Envelope, 1.2x. El Envelope es
// un gancho opcional para overlays de scoring o cambios de material de
// fondo sin tocar el World (ver documento de diseño, sección 2.1).
constexpr G4double kWorldScale = 1.4;
constexpr G4double kEnvelopeScale = 1.2;

// Defaults de la geometría
constexpr G4double kDefaultThickness = 5.0 * CLHEP::mm;
constexpr G4double kDefaultSizeXY = 10.0 * CLHEP::cm;
constexpr G4int kDefaultNumberOfLayers = 1;
const G4String kDefaultMaterial = "G4_WATER";

// Chequeo de overlaps en todos los G4PVPlacement (requisito de diseño).
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
  // Al reconstruir la geometría tras un cambio por macro hay que limpiar los
  // stores; si no, los volúmenes viejos quedan colgando y Geant4 aborta.
  // Antes de limpiar, desregistramos el volumen raíz de la región para no
  // dejar un puntero colgante dentro de "AbsorberRegion".
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
  // World: caja de aire 1.4x el slab. Nota física: el protón primario nace en
  // aire a ~0.15*espesor de la cara de entrada (ver PrimaryGeneratorAction);
  // la pérdida de energía en ese trayecto de aire es despreciable (<1 keV
  // para las energías de este estudio), pero si se quisiera eliminar por
  // completo bastaría cambiar el material del World a G4_Galactic aquí.
  // --------------------------------------------------------------------------
  const G4double worldHalfXY = 0.5 * kWorldScale * fSizeXY;
  const G4double worldHalfZ = 0.5 * kWorldScale * fThickness;
  auto* worldSolid = new G4Box("World", worldHalfXY, worldHalfXY, worldHalfZ);
  fWorldLV = new G4LogicalVolume(worldSolid, air, "World");
  fWorldPV = new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLV, "World", nullptr,
                               false, 0, kCheckOverlaps);

  // --------------------------------------------------------------------------
  // Envelope: mismo material que el World, 1.2x el slab.
  // --------------------------------------------------------------------------
  const G4double envHalfXY = 0.5 * kEnvelopeScale * fSizeXY;
  const G4double envHalfZ = 0.5 * kEnvelopeScale * fThickness;
  auto* envSolid = new G4Box("Envelope", envHalfXY, envHalfXY, envHalfZ);
  fEnvelopeLV = new G4LogicalVolume(envSolid, air, "Envelope");
  new G4PVPlacement(nullptr, G4ThreeVector(), fEnvelopeLV, "Envelope", fWorldLV,
                    false, 0, kCheckOverlaps);

  // --------------------------------------------------------------------------
  // Slab absorbedor: G4Box con el eje del haz en z.
  // Si fNumberOfLayers > 1, el slab contenedor se subdivide en capas iguales
  // mediante G4PVReplica a lo largo de kZAxis; cada capa recibe un índice de
  // copia (copyNo) que SteppingAction usa para el scoring por profundidad.
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
    // La capa debe tener el mismo material que el contenedor (la réplica
    // llena completamente a su madre).
    fLayerLV = new G4LogicalVolume(layerSolid, fAbsorberMaterial, "Layer");
    // G4PVReplica no acepta bandera de chequeo de overlaps: la subdivisión
    // regular es exacta por construcción.
    new G4PVReplica("Layer", fLayerLV, fAbsorberLV, kZAxis, fNumberOfLayers,
                    layerThickness);
  } else {
    fLayerLV = nullptr;
  }

  // --------------------------------------------------------------------------
  // Región propia del slab, para cortes de producción locales
  // (/physics/absorberCut, ver PhysicsList::SetCuts()).
  // --------------------------------------------------------------------------
  fRegion = G4RegionStore::GetInstance()->GetRegion(kRegionName, false);
  if (fRegion == nullptr) {
    fRegion = new G4Region(kRegionName);
  }
  fRegion->AddRootLogicalVolume(fAbsorberLV);

  // --- Atributos de visualización mínimos ---
  fWorldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* envVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.1));
  envVis->SetForceWireframe(true);
  fEnvelopeLV->SetVisAttributes(envVis);
  auto* slabVis = new G4VisAttributes(G4Colour(0.2, 0.4, 0.8, 0.3));
  slabVis->SetForceSolid(true);
  fAbsorberLV->SetVisAttributes(slabVis);

  G4cout << "### DetectorConstruction: slab de "
         << fAbsorberMaterial->GetName() << ", espesor "
         << G4BestUnit(fThickness, "Length") << ", seccion "
         << G4BestUnit(fSizeXY, "Length") << " x "
         << G4BestUnit(fSizeXY, "Length") << ", " << fNumberOfLayers
         << " capa(s)" << G4endl;

  return fWorldPV;
}

// ----------------------------------------------------------------------------
// Setters de mensajería. La reconstrucción efectiva la dispara
// DetectorMessenger tras llamar a estos métodos (ver DetectorMessenger.cc).
// ----------------------------------------------------------------------------
void DetectorConstruction::SetAbsorberMaterial(const G4String& nistName)
{
  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(nistName);
  if (mat == nullptr) {
    G4cerr << "### DetectorConstruction: material NIST desconocido '"
           << nistName << "'. Se conserva "
           << fAbsorberMaterial->GetName() << G4endl;
    return;
  }
  fAbsorberMaterial = mat;
}

void DetectorConstruction::SetThickness(G4double value)
{
  if (value <= 0.) {
    G4cerr << "### DetectorConstruction: espesor invalido, se ignora" << G4endl;
    return;
  }
  fThickness = value;
}

void DetectorConstruction::SetSizeXY(G4double value)
{
  if (value <= 0.) {
    G4cerr << "### DetectorConstruction: sizeXY invalido, se ignora" << G4endl;
    return;
  }
  fSizeXY = value;
}

void DetectorConstruction::SetNumberOfLayers(G4int n)
{
  if (n < 1) {
    G4cerr << "### DetectorConstruction: numberOfLayers < 1, se ignora"
           << G4endl;
    return;
  }
  fNumberOfLayers = n;
}
