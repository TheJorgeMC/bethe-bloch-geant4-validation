// ============================================================================
// DetectorConstruction.hh
// Geometry: World (G4_AIR) -> Envelope (G4_AIR) -> absorber slab (G4Box),
// with thickness, material, transverse size and number of layers configurable
// via UI commands (/absorber/...). The slab defines the G4Region
// "AbsorberRegion" so that region-specific production cuts can be applied
// (see PhysicsList).
// ============================================================================
#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  // --- Setters invoked by DetectorMessenger (/absorber/... commands) ---
  void SetAbsorberMaterial(const G4String& nistName);
  void SetThickness(G4double value);
  void SetSizeXY(G4double value);
  void SetNumberOfLayers(G4int n);

  // --- Getters used by PrimaryGeneratorAction, RunAction and Stepping ---
  G4Material* GetAbsorberMaterial() const { return fAbsorberMaterial; }
  G4double GetThickness() const { return fThickness; }
  G4double GetSizeXY() const { return fSizeXY; }
  G4int GetNumberOfLayers() const { return fNumberOfLayers; }

  // Returns true if the logical volume belongs to the absorber slab (the
  // container or one of its replicated layers). Used by SteppingAction to
  // decide whether a step occurs inside the volume of interest.
  G4bool IsAbsorber(const G4LogicalVolume* lv) const
  {
    return (lv != nullptr) && (lv == fAbsorberLV || lv == fLayerLV);
  }

  static constexpr const char* kRegionName = "AbsorberRegion";

 private:
  void CleanupGeometry();

  // --- Configurable parameters (defaults set in the constructor) ---
  G4Material* fAbsorberMaterial = nullptr;
  G4double fThickness;       // total slab thickness along z
  G4double fSizeXY;          // side of the (square) transverse cross section
  G4int fNumberOfLayers;     // 1 = single slab; N>1 = replicated along z

  // --- Volumes and region ---
  G4LogicalVolume* fWorldLV = nullptr;
  G4VPhysicalVolume* fWorldPV = nullptr;
  G4LogicalVolume* fEnvelopeLV = nullptr;
  G4LogicalVolume* fAbsorberLV = nullptr;  // slab container
  G4LogicalVolume* fLayerLV = nullptr;     // replicated layer (only if N>1)
  G4Region* fRegion = nullptr;

  DetectorMessenger* fMessenger = nullptr;
};

#endif
