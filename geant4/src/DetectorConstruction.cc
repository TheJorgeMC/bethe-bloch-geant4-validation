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
#include "G4RotationMatrix.hh"
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
#include <cmath>

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

// Slab colour -- darkened per request (was (0.2, 0.4, 0.8, 0.3)): same hue,
// noticeably deeper blue so it reads clearly against the new white
// background set in vis.mac (/vis/viewer/set/background white). Used for
// both fAbsorberLV (numberOfLayers == 1) and fLayerLV (numberOfLayers > 1,
// see the vis-attributes block in Construct()).
const G4Colour kSlabColour(0.06, 0.18, 0.45, 0.55);

// ----------------------------------------------------------------------------
// Decorative logos (LIFI / CUCEI banner / CUCEI "R") — visualization only.
//
// WHY THIS IS SAFE FOR THE PHYSICS VALIDATION:
//  - Gated by fShowLogos (default false, set via /detector/showLogos — see
//    DetectorMessenger). Every production/energy-scan run that never issues
//    that command gets EXACTLY the same World/Envelope/Absorber geometry as
//    before this change, byte-for-byte: the World is only enlarged, and the
//    logo volumes only built, when fShowLogos is true.
//  - The logos live inside a dedicated, invisible "LogoBillboard" container
//    that is itself a child of the WORLD logical volume, a SIBLING of
//    Envelope — not a descendant of Absorber — so nothing here is ever added
//    to fRegion ("AbsorberRegion": see fRegion->AddRootLogicalVolume below,
//    unchanged) and therefore cannot affect /physics/absorberCut or any
//    other region-scoped production cut.
//  - Each logo shape gets its own freshly-constructed G4LogicalVolume (name
//    "<group>_<index>_LV" from LogoGeometryBuilder.cc) that is never passed
//    to SetSensitiveDetector anywhere in this file. Geant4 scoring is
//    opt-in per logical volume (not by geometric proximity), so these
//    volumes cannot produce hits in any G4MultiFunctionalDetector this
//    project already has configured, regardless of where they are placed.
//
// LAYOUT: side by side (NOT stacked in z), left to right in the order
// requested: LIFI, CUCEI banner, CUCEI "R". Scaled to 40% of their native
// size (kLogoScale) via G4ExtrudedSolid's own scale1/scale2 parameters (see
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

// Centre-x of each logo, left to right, with the whole row centred on its
// OWN local x=0 (see kLogoFaceDirection below: the row is built flat, in
// its own local x-y plane, then rotated+translated as a single rigid
// "billboard" -- the per-logo positions here never change) -- LIFI first,
// then the banner, then the "R", exactly the order requested.
constexpr G4double kLogoLIFICenterXMm =
    -0.5 * kLogoRowWidthMm + 0.5 * kLogoLIFIWidthMm;                    // -66.0
constexpr G4double kLogoBannerCenterXMm =
    kLogoLIFICenterXMm + 0.5 * kLogoLIFIWidthMm + kLogoGapMm +
    0.5 * kLogoBannerWidthMm;                                           // +16.0
constexpr G4double kLogoCUCEIRCenterXMm =
    kLogoBannerCenterXMm + 0.5 * kLogoBannerWidthMm + kLogoGapMm +
    0.5 * kLogoCUCEIRWidthMm;                                           // +82.0

constexpr G4double kLogoMaxHeightMm =
    std::max({kLogoLIFIHeightMm, kLogoBannerHeightMm, kLogoCUCEIRHeightMm});

// Fixed (NOT thickness-scaled) clearance past the Envelope's outer face,
// and fixed margin added on top of the logos' own footprint when sizing
// the World/billboard container. Both independent of the current
// (messenger-settable) slab thickness, so the logos never touch the
// Envelope no matter what /absorber/thickness is set to.
constexpr G4double kLogoClearanceMm = 3.0;
constexpr G4double kLogoWorldMarginMm = 5.0;

// Fixed vertical position of the logo row: -40% of the DEFAULT (not the
// current messenger-settable) sizeXY, so the row sits near the bottom of
// the visualization scene at a position that does not silently shift if
// someone runs /absorber/sizeXY before /detector/showLogos. Requested
// value: -0.4 * kDefaultSizeXY.
constexpr G4double kLogoBottomYMm = -0.4 * kDefaultSizeXY;

// Extra empty space (vis-only, see the fShowLogos block in Construct())
// so the incoming proton track has visible room to travel before reaching
// the entrance face. NOTE this only creates the *space* -- see the caveat
// where it's used below.
constexpr G4double kVisBeamApproachMm = 150.0;

// Camera direction: "a vector between the positive x, y, z axes", i.e.
// equidistant from all three -- (1,1,1) normalized. Matches
// /vis/viewer/set/viewpointVector 1 1 1 in vis.mac. The logo billboard is
// rotated so its own local z-axis (the flat 2D face's normal, the axis
// G4ExtrudedSolid extrudes along) points along this direction, so the row
// of logos appears fully face-on when viewed from that viewpoint.
const G4ThreeVector kLogoFaceDirection = G4ThreeVector(0.759055, 0.115654, 0.640671).unit();

// Angle/axis of the fixed billboard tilt, derived once and reused both for
// the G4PVPlacement rotation and for the corner-extent math below. Geant4
// convention: the rotation matrix passed to G4PVPlacement maps MOTHER-frame
// vectors into the daughter's LOCAL frame (local = pRot * mother) -- NOT
// the intuitive "rotate the object by this matrix".
const G4ThreeVector kBillboardAxis =
    G4ThreeVector(0., 0., 1.).cross(kLogoFaceDirection).unit();
const G4double kBillboardAngle =
    std::acos(G4ThreeVector(0., 0., 1.).dot(kLogoFaceDirection));

// FORWARD rotation: local -> world (local z ends up pointing along
// kLogoFaceDirection). This is the physical rotation the billboard
// actually undergoes -- the inverse/transpose of what G4PVPlacement wants
// (see MakeLogoBillboardRotation below). Used only for the corner-extent
// math in ComputeBillboardExtents(), never handed to G4PVPlacement itself.
G4RotationMatrix ForwardBillboardRotation()
{
  G4RotationMatrix rot;
  rot.rotate(kBillboardAngle, kBillboardAxis);
  return rot;
}

// Builds the rotation matrix to hand to G4PVPlacement for the billboard
// (the inverse of ForwardBillboardRotation() -- same axis, negated angle).
G4RotationMatrix* MakeLogoBillboardRotation()
{
  auto* rot = new G4RotationMatrix();
  rot->rotate(-kBillboardAngle, kBillboardAxis);
  return rot;
}

// World-frame extents of the billboard box's 8 corners AFTER
// ForwardBillboardRotation() is applied but BEFORE any translation --
// i.e. "how far does the tilted panel reach in each direction, measured
// from its own (not yet placed) centre". Used by Construct() to work out
// the TIGHT minimum z-offset that clears the Envelope (rather than a
// generic, much larger bounding-sphere guess) and the World half-extents
// actually needed to contain the tilted billboard.
//
// Why checking minZ alone is enough to guarantee no overlap: the box is
// convex, so if the offset is chosen such that its closest corner (minZ)
// sits at world z = envHalfZ + clearance, EVERY other corner has a z
// greater than that too -- the whole box then lies entirely in the
// half-space z > envHalfZ, which cannot intersect the Envelope (confined
// to |z| <= envHalfZ), regardless of the box's x/y footprint at that z.
struct BillboardExtents
{
  G4double minZ = 0.;
  G4double maxZ = 0.;
  G4double maxAbsX = 0.;
  G4double maxAbsY = 0.;
};

BillboardExtents ComputeBillboardExtents(G4double halfX, G4double halfY,
                                          G4double halfZ)
{
  const G4RotationMatrix rot = ForwardBillboardRotation();
  BillboardExtents ext;
  G4bool first = true;
  for (G4double sx : {-1., 1.}) {
    for (G4double sy : {-1., 1.}) {
      for (G4double sz : {-1., 1.}) {
        const G4ThreeVector corner =
            rot * G4ThreeVector(sx * halfX, sy * halfY, sz * halfZ);
        if (first) {
          ext.minZ = ext.maxZ = corner.z();
          first = false;
        } else {
          ext.minZ = std::min(ext.minZ, corner.z());
          ext.maxZ = std::max(ext.maxZ, corner.z());
        }
        ext.maxAbsX = std::max(ext.maxAbsX, std::abs(corner.x()));
        ext.maxAbsY = std::max(ext.maxAbsY, std::abs(corner.y()));
      }
    }
  }
  return ext;
}

// Half-extents of the invisible billboard container -- shared by
// PlaceLogos() (to actually build the box) and by Construct() (to run
// ComputeBillboardExtents() before the box exists, since it only needs
// numbers, not the G4Box object itself).
G4double BillboardHalfX() { return 0.5 * kLogoRowWidthMm + kLogoWorldMarginMm; }
G4double BillboardHalfY() { return 0.5 * kLogoMaxHeightMm + kLogoWorldMarginMm; }
G4double BillboardHalfZ() { return kLogoHalfThicknessMm + kLogoWorldMarginMm; }

// Places the three logos side by side inside their own local x-y plane
// (z=0 for all of them, laid out along local x only), left to right: LIFI,
// CUCEI banner, CUCEI "R" — as requested. That flat row is wrapped in an
// invisible "LogoBillboard" container which is placed into 'worldLV' at
// 'billboardCentre' with the rotation from MakeLogoBillboardRotation(), so
// the row is tilted as one rigid unit to face the camera. 'billboardCentre'
// is already validated by the caller (Construct(), via
// ComputeBillboardExtents()) to keep the whole tilted billboard clear of
// the Envelope.
void PlaceLogos(G4LogicalVolume* worldLV, const G4ThreeVector& billboardCentre,
                 G4bool checkOverlaps)
{
  auto* nist = G4NistManager::Instance();
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  // Purely decorative: any solid, visually opaque material works — it does
  // not need to match a "real" detector material, and (see the note above)
  // cannot influence dose/dE/dx scoring regardless of what it's made of.
  G4Material* logoMaterial = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

  const G4double halfX = BillboardHalfX();
  const G4double halfY = BillboardHalfY();
  const G4double halfZ = BillboardHalfZ();
  auto* billboardSolid =
      new G4Box("LogoBillboard", halfX * mm, halfY * mm, halfZ * mm);
  auto* billboardLV = new G4LogicalVolume(billboardSolid, air, "LogoBillboard");
  billboardLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4RotationMatrix* billboardRot = MakeLogoBillboardRotation();
  new G4PVPlacement(billboardRot, billboardCentre, billboardLV, "LogoBillboard",
                     worldLV, false, 0, checkOverlaps);

  BuildLogoGeometry(BuildLIFILogoData(), kLogoHalfThicknessMm, logoMaterial,
                    billboardLV,
                    G4ThreeVector(kLogoLIFICenterXMm * mm, 0., 0.),
                    checkOverlaps, /*copyNoStart=*/0, kLogoScale);
  BuildLogoGeometry(BuildCUCEIBannerLogoData(), kLogoHalfThicknessMm,
                    logoMaterial, billboardLV,
                    G4ThreeVector(kLogoBannerCenterXMm * mm, 0., 0.),
                    checkOverlaps, /*copyNoStart=*/0, kLogoScale);
  BuildLogoGeometry(BuildCUCEI_R_LogoData(), kLogoHalfThicknessMm,
                    logoMaterial, billboardLV,
                    G4ThreeVector(kLogoCUCEIRCenterXMm * mm, 0., 0.),
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
  // logos-disabled run gets the exact same World it always did. The
  // billboard's centre is computed here (before the Envelope solid exists)
  // since envHalfXY/envHalfZ are pure numbers, not objects.
  G4ThreeVector logoBillboardCentre;
  if (fShowLogos) {
    const G4double envHalfXYEst = 0.5 * kEnvelopeScale * fSizeXY;
    const G4double envHalfZEst = 0.5 * kEnvelopeScale * fThickness;

    // TIGHT (not a generic bounding-sphere guess) placement: work out
    // exactly how far the tilted billboard's corners reach, in world axes,
    // around its own (not yet placed) centre -- see the doc comment on
    // ComputeBillboardExtents() above for why checking minZ/maxZ this way
    // is sufficient to guarantee no overlap with the Envelope, regardless
    // of the billboard's x/y footprint.
    const BillboardExtents ext = ComputeBillboardExtents(
        BillboardHalfX(), BillboardHalfY(), BillboardHalfZ());

    const G4double logoCentreXMm = 0.;
    const G4double logoCentreYMm = kLogoBottomYMm;
    // Push the centre along +z just enough that the NEAREST corner
    // (ext.minZ, always <= 0) lands exactly 'kLogoClearanceMm' past the
    // Envelope's +z face.
    const G4double logoCentreZMm = envHalfZEst + kLogoClearanceMm - ext.minZ;
    logoBillboardCentre = G4ThreeVector(logoCentreXMm * mm, logoCentreYMm * mm,
                                        logoCentreZMm * mm);

    const G4double neededHalfXY =
        std::max(std::abs(logoCentreXMm) + ext.maxAbsX,
                 std::abs(logoCentreYMm) + ext.maxAbsY) + kLogoWorldMarginMm;
    const G4double neededHalfZ =
        logoCentreZMm + ext.maxZ + kLogoWorldMarginMm;

    worldHalfXY = std::max(worldHalfXY, neededHalfXY * mm);
    worldHalfZ = std::max(worldHalfZ, neededHalfZ * mm);

    // Extra room so the incoming proton track is visibly long before it
    // reaches the entrance face. The World is symmetric about the origin,
    // so growing worldHalfZ opens up space on BOTH the upstream (-z, where
    // the gun fires from) and downstream (+z) sides.
    // CAVEAT: this only creates the *space*. If PrimaryGeneratorAction
    // places the gun at a position derived from the World's own extent,
    // the beam will now visibly start further out. If instead it uses a
    // FIXED offset from the entrance face (as noted above,
    // "~0.15*thickness upstream" — independent of World size), growing the
    // World alone will not move the gun; PrimaryGeneratorAction.cc would
    // need a matching change too.
    worldHalfZ = std::max(worldHalfZ, kVisBeamApproachMm * mm);
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
  // chance to render in that case. Give the Layer LV the SAME darker blue,
  // plus SetForceAuxEdgeVisible(true) so the boundary between adjacent
  // layers is explicitly drawn instead of the whole stack merging into one
  // solid blue block — this is what makes the per-layer division visible.
  //
  // IMPORTANT LIMITATION (not fixable from this file): G4VisAttributes has
  // a single RGBA colour, used for both the filled surface (ForceSolid) AND
  // any wireframe/auxiliary edges drawn on it — there is no separate
  // "edge colour" API. The colour of the layer-boundary LINES themselves is
  // controlled by the VIEWER's global default line colour, i.e. the
  // /vis/viewer/set/defaultLineColour command in vis.mac, not by anything a
  // G4LogicalVolume's vis attributes can override. Your current vis.mac
  // sets that to black; to get dark-grey layer lines (clearly distinct from
  // the darker blue fill against the new white background), change that one
  // line in vis.mac to e.g.:
  //   /vis/viewer/set/defaultLineColour grey
  // or an explicit dark grey: /vis/viewer/set/defaultLineColour 0.3 0.3 0.3
  if (fLayerLV != nullptr) {
    auto* layerVis = new G4VisAttributes(kSlabColour);
    layerVis->SetForceSolid(true);
    layerVis->SetForceAuxEdgeVisible(true);
    fLayerLV->SetVisAttributes(layerVis);
  }

  // --------------------------------------------------------------------------
  // Decorative logos (LIFI / CUCEI banner / CUCEI "R") — visualization only.
  // See the PlaceLogos()/kLogo* block above (anonymous namespace) for the
  // full rationale. Placed LAST, as children of fWorldLV (via the invisible
  // "LogoBillboard" container), strictly after the physics geometry
  // (World/Envelope/Absorber/Layer/Region) is fully built, so their
  // presence can never shadow or reorder anything above.
  // --------------------------------------------------------------------------
  if (fShowLogos) {
    PlaceLogos(fWorldLV, logoBillboardCentre, kCheckOverlaps);
    G4cout << "### DetectorConstruction: decorative logos placed as a "
              "camera-facing billboard centred at "
           << G4BestUnit(logoBillboardCentre, "Length") << " (scale "
           << (kLogoScale * 100.) << "%, left-to-right: LIFI, CUCEI banner, "
              "CUCEI R, facing (1,1,1); World half-extents: xy = "
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
  // Purely decorative flag (see the PlaceLogos()/kLogo* block above for why
  // this cannot affect physics): toggling it and re-running the geometry
  // (via ReinitializeGeometry(), triggered by DetectorMessenger) only adds
  // or removes the logo billboard and, if needed, enlarges the World.
  fShowLogos = value;
}