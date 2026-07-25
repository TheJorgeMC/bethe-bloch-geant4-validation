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

// --- Decorative logo geometry (visualization only) -------------------------
// See LogoIntegration section below for why this lives here rather than in
// its own translation unit: it only needs four headers and ~30 lines, and
// keeping it next to Construct() makes the "logos never touch the physics
// geometry" guarantee easy to audit in one file.
#include "LIFILogoData.hh"
#include "CUCEI_R_LogoData.hh"
#include "CUCEIBannerLogoData.hh"
#include "LogoGeometryBuilder.hh"

#include <algorithm>

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

// Slab colour (same blue used whether the slab is a single block or
// segmented into layers -- see the vis-attributes block in Construct()).
const G4Colour kSlabColour(0.2, 0.4, 0.8, 0.3);

// ----------------------------------------------------------------------------
// Decorative logos (LIFI / CUCEI banner / CUCEI "R") — visualization only.
//
// WHY THIS IS SAFE FOR THE PHYSICS VALIDATION:
//  - Gated by fShowLogos (default false, set via /detector/showLogos — see
//    DetectorMessenger). Every production/energy-scan run that never issues
//    that command gets EXACTLY the same World/Envelope/Absorber geometry as
//    before this change, byte-for-byte: the World is only enlarged, and the
//    logo volumes only built, when fShowLogos is true.
//  - The logos are children of the WORLD logical volume, a SIBLING of
//    Envelope — not descendants of Absorber — so they are NOT added to
//    fRegion ("AbsorberRegion": see fRegion->AddRootLogicalVolume(fAbsorberLV)
//    below, unchanged) and therefore cannot affect /physics/absorberCut or
//    any other region-scoped production cut.
//  - Each logo shape gets its own freshly-constructed G4LogicalVolume (name
//    "<group>_<index>_LV" from LogoGeometryBuilder.cc) that is never passed
//    to SetSensitiveDetector anywhere in this file. Geant4 scoring is
//    opt-in per logical volume (not by geometric proximity), so these
//    volumes cannot produce hits in any G4MultiFunctionalDetector this
//    project already has configured, regardless of where they are placed.
//
// LAYOUT: side by side (NOT stacked in z), all at the SAME z-plane just
// downstream of the Envelope, left to right in the order requested:
// LIFI, CUCEI banner, CUCEI "R". Scaled to 40% of their native size
// (kLogoScale) via G4ExtrudedSolid's own scale1/scale2 parameters (see
// LogoGeometryBuilder.cc) -- the vertex data itself is untouched.
//
// Native (unscaled) bounding footprints, from *LogoData.hh (local origin =
// each logo's own image centre): LIFI 160.00 x 50.46 mm, CUCEI banner
// 200.00 x 80.00 mm, CUCEI "R" 80.00 x 113.15 mm.
constexpr G4double kLogoScale = 0.4;             // 40% of native size
constexpr G4double kLogoHalfThicknessMm = 1.0;   // 2 mm-thick plaques
                                                  // (thickness is NOT scaled)

constexpr G4double kLogoLIFIWidthMm = 160.00 * kLogoScale;    // 64.00
constexpr G4double kLogoLIFIHeightMm = 50.46 * kLogoScale;    // 20.184
constexpr G4double kLogoBannerWidthMm = 200.00 * kLogoScale;  // 80.00
constexpr G4double kLogoBannerHeightMm = 80.00 * kLogoScale;  // 32.00
constexpr G4double kLogoCUCEIRWidthMm = 80.00 * kLogoScale;   // 32.00
constexpr G4double kLogoCUCEIRHeightMm = 113.15 * kLogoScale; // 45.26

// Gap between adjacent logo bounding boxes in the row, and the row's total
// width (LIFI + gap + Banner + gap + CUCEI R).
constexpr G4double kLogoGapMm = 10.0;
constexpr G4double kLogoRowWidthMm = kLogoLIFIWidthMm + kLogoGapMm +
                                     kLogoBannerWidthMm + kLogoGapMm +
                                     kLogoCUCEIRWidthMm;  // 196.0

// Centre-x of each logo, left to right, with the whole row centred on the
// beam axis (x=0) -- LIFI first, then the banner, then the "R", exactly the
// order requested.
constexpr G4double kLogoLIFICenterXMm =
    -0.5 * kLogoRowWidthMm + 0.5 * kLogoLIFIWidthMm;                    // -66.0
constexpr G4double kLogoBannerCenterXMm =
    kLogoLIFICenterXMm + 0.5 * kLogoLIFIWidthMm + kLogoGapMm +
    0.5 * kLogoBannerWidthMm;                                           // +16.0
constexpr G4double kLogoCUCEIRCenterXMm =
    kLogoBannerCenterXMm + 0.5 * kLogoBannerWidthMm + kLogoGapMm +
    0.5 * kLogoCUCEIRWidthMm;                                           // +82.0

// Fixed (NOT thickness-scaled) clearance past the Envelope's outer face,
// and fixed margin added on top of the logos' own footprint when sizing
// the World. Both independent of the current (messenger-settable) slab
// thickness, so the logos never touch the Envelope no matter what
// /absorber/thickness is set to.
constexpr G4double kLogoClearanceMm = 10.0;
constexpr G4double kLogoWorldMarginMm = 5.0;

// Places the three logos side by side (single z-plane, different x), left
// to right: LIFI, CUCEI banner, CUCEI "R" — as requested. z already
// computed and validated to fit inside the World by the caller.
void PlaceLogos(G4LogicalVolume* worldLV, G4double zLogo, G4bool checkOverlaps)
{
  // Purely decorative: any solid, visually opaque material works — it does
  // not need to match a "real" detector material, and (see the note above)
  // cannot influence dose/dE/dx scoring regardless of what it's made of.
  G4Material* logoMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");

  BuildLogoGeometry(BuildLIFILogoData(), kLogoHalfThicknessMm, logoMaterial,
                    worldLV,
                    G4ThreeVector(kLogoLIFICenterXMm * mm, 0., zLogo),
                    checkOverlaps, /*copyNoStart=*/0, kLogoScale);
  BuildLogoGeometry(BuildCUCEIBannerLogoData(), kLogoHalfThicknessMm,
                    logoMaterial, worldLV,
                    G4ThreeVector(kLogoBannerCenterXMm * mm, 0., zLogo),
                    checkOverlaps, /*copyNoStart=*/0, kLogoScale);
  BuildLogoGeometry(BuildCUCEI_R_LogoData(), kLogoHalfThicknessMm,
                    logoMaterial, worldLV,
                    G4ThreeVector(kLogoCUCEIRCenterXMm * mm, 0., zLogo),
                    checkOverlaps, /*copyNoStart=*/0, kLogoScale);
}
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
  // World: air box, 1.4x the slab (or larger still if the decorative logos
  // are enabled — see below). Physics note: the primary proton is born in
  // air at ~0.15*thickness upstream of the entrance face (see
  // PrimaryGeneratorAction); the energy loss along that air path is
  // negligible (<1 keV for the energies of this study). To remove it
  // entirely, change the World material to G4_Galactic here.
  // --------------------------------------------------------------------------
  G4double worldHalfXY = 0.5 * kWorldScale * fSizeXY;
  G4double worldHalfZ = 0.5 * kWorldScale * fThickness;

  // If logos are enabled, grow the World ONLY as much as needed to contain
  // them, via max(), not by replacing the physics-derived size: a
  // logos-disabled run gets the exact same World it always did. The z
  // position is computed here (before the Envelope solid exists) since
  // envHalfZ is a pure number, not an object.
  G4double logoZ = 0.;
  if (fShowLogos) {
    const G4double envHalfZ = 0.5 * kEnvelopeScale * fThickness;
    logoZ = envHalfZ + (kLogoClearanceMm + kLogoHalfThicknessMm) * mm;
    const G4double neededHalfZ =
        logoZ + (kLogoHalfThicknessMm + kLogoWorldMarginMm) * mm;
    // Lateral half-extent needed: the wider of (a) half the total row
    // width (x direction) and (b) the tallest single logo's half-height
    // (y direction, since every logo is centred at y=0).
    const G4double neededHalfXY =
        std::max(0.5 * kLogoRowWidthMm,
                 0.5 * std::max({kLogoLIFIHeightMm, kLogoBannerHeightMm,
                                 kLogoCUCEIRHeightMm})) +
        kLogoWorldMarginMm;
    worldHalfXY = std::max(worldHalfXY, neededHalfXY * mm);
    worldHalfZ = std::max(worldHalfZ, neededHalfZ);
  }

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

  auto* slabVis = new G4VisAttributes(kSlabColour);
  slabVis->SetForceSolid(true);
  fAbsorberLV->SetVisAttributes(slabVis);

  // When segmented (fNumberOfLayers > 1), it's the Layer REPLICA that is
  // actually drawn (the replicas fully cover the Absorber container it
  // lives in), so the Absorber's own vis attributes above never get a
  // chance to render in that case. Give the Layer LV the SAME blue, plus
  // SetForceAuxEdgeVisible(true) so the boundary between adjacent layers
  // is explicitly drawn instead of the whole stack merging into one solid
  // blue block — this is what makes the per-layer division visible again.
  if (fLayerLV != nullptr) {
    auto* layerVis = new G4VisAttributes(kSlabColour);
    layerVis->SetForceSolid(true);
    layerVis->SetForceAuxEdgeVisible(true);
    fLayerLV->SetVisAttributes(layerVis);
  }

  // --------------------------------------------------------------------------
  // Decorative logos (LIFI / CUCEI banner / CUCEI "R") — visualization only.
  // See the PlaceLogos()/kLogo* block above (anonymous namespace) for the
  // full rationale. Placed LAST, as children of fWorldLV, strictly after
  // the physics geometry (World/Envelope/Absorber/Layer/Region) is fully
  // built, so their presence can never shadow or reorder anything above.
  // --------------------------------------------------------------------------
  if (fShowLogos) {
    PlaceLogos(fWorldLV, logoZ, kCheckOverlaps);
    G4cout << "### DetectorConstruction: decorative logos placed at z = "
           << G4BestUnit(logoZ, "Length") << " (scale "
           << (kLogoScale * 100.) << "%, left-to-right: LIFI, CUCEI banner, "
              "CUCEI R; World half-extents: xy = "
           << G4BestUnit(worldHalfXY, "Length")
           << ", z = " << G4BestUnit(worldHalfZ, "Length") << ")" << G4endl;
  }

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

void DetectorConstruction::SetShowLogos(G4bool value)
{
  // No validation needed (unlike thickness/sizeXY/numberOfLayers): both
  // true and false are always valid. Purely decorative — see the
  // PlaceLogos()/kLogo* block above for why this cannot affect physics
  // results even though it changes the World's size when true.
  fShowLogos = value;
}