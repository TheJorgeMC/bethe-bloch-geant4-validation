// ============================================================================
// DetectorConstruction.hh
// Geometría: World (G4_AIR) -> Envelope (G4_AIR) -> Slab absorbedor (G4Box),
// con espesor, material, sección transversal y número de capas configurables
// por mensajería (/absorber/...). El slab define la G4Region "AbsorberRegion"
// para poder aplicar cortes de producción específicos (ver PhysicsList).
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

  // --- Setters invocados por DetectorMessenger (comandos /absorber/...) ---
  void SetAbsorberMaterial(const G4String& nistName);
  void SetThickness(G4double value);
  void SetSizeXY(G4double value);
  void SetNumberOfLayers(G4int n);

  // --- Getters usados por PrimaryGeneratorAction, RunAction y Stepping ---
  G4Material* GetAbsorberMaterial() const { return fAbsorberMaterial; }
  G4double GetThickness() const { return fThickness; }
  G4double GetSizeXY() const { return fSizeXY; }
  G4int GetNumberOfLayers() const { return fNumberOfLayers; }

  // Devuelve true si el volumen lógico pertenece al slab absorbedor
  // (el contenedor o una de sus capas replicadas). Lo usa SteppingAction
  // para decidir si un paso ocurre dentro del volumen de interés.
  G4bool IsAbsorber(const G4LogicalVolume* lv) const
  {
    return (lv != nullptr) && (lv == fAbsorberLV || lv == fLayerLV);
  }

  static constexpr const char* kRegionName = "AbsorberRegion";

 private:
  void CleanupGeometry();

  // --- Parámetros configurables (defaults fijados en el constructor) ---
  G4Material* fAbsorberMaterial = nullptr;
  G4double fThickness;       // espesor total del slab a lo largo de z
  G4double fSizeXY;          // lado de la sección transversal (cuadrada)
  G4int fNumberOfLayers;     // 1 = slab único; N>1 = réplica en z

  // --- Volúmenes y región ---
  G4LogicalVolume* fWorldLV = nullptr;
  G4VPhysicalVolume* fWorldPV = nullptr;
  G4LogicalVolume* fEnvelopeLV = nullptr;
  G4LogicalVolume* fAbsorberLV = nullptr;  // contenedor del slab
  G4LogicalVolume* fLayerLV = nullptr;     // capa replicada (solo si N>1)
  G4Region* fRegion = nullptr;

  DetectorMessenger* fMessenger = nullptr;
};

#endif
