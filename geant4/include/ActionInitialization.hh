// ============================================================================
// ActionInitialization.hh
// Standard wiring of the user actions for sequential and multi-threaded mode.
// ============================================================================
#ifndef ActionInitialization_hh
#define ActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization
{
 public:
  explicit ActionInitialization(const DetectorConstruction* det);
  ~ActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;

 private:
  const DetectorConstruction* fDetector;
};

#endif
