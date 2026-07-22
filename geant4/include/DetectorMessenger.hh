// ============================================================================
// DetectorMessenger.hh
// UI commands for the slab geometry:
//   /absorber/material <NIST name>
//   /absorber/thickness <value> <unit>
//   /absorber/sizeXY <value> <unit>
//   /absorber/numberOfLayers <N>
// ============================================================================
#ifndef DetectorMessenger_hh
#define DetectorMessenger_hh 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class DetectorMessenger : public G4UImessenger
{
 public:
  explicit DetectorMessenger(DetectorConstruction* det);
  ~DetectorMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String newValue) override;

 private:
  DetectorConstruction* fDetector;

  G4UIdirectory* fDirectory = nullptr;
  G4UIcmdWithAString* fMaterialCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fThicknessCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fSizeXYCmd = nullptr;
  G4UIcmdWithAnInteger* fLayersCmd = nullptr;
};

#endif
