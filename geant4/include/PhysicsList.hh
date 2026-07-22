// ============================================================================
// PhysicsList.hh
// Minimal modular list: ONLY standard electromagnetic physics
// (G4EmStandardPhysics_option4 by default; option3 selectable as a
// sensitivity check via command-line argument in slab.cc:
//   ./slab run.mac 3
// ).
//
// Deliberately WITHOUT hadronic physics: for protons crossing a thin slab,
// inelastic nuclear reactions are a second-order effect for the dE/dx
// validation and would only add a primary-loss channel and computing time
// (decision documented in the design, section 4).
//
// Slab-specific production cut:
//   /physics/absorberCut <value> <unit>   (default 1 mm)
// applied to the G4Region "AbsorberRegion" in SetCuts(). This cut sets the
// threshold T_cut for explicit delta-ray production and therefore defines
// which restricted dE/dx the continuous ionisation process uses (design,
// section 5).
// ============================================================================
#ifndef PhysicsList_hh
#define PhysicsList_hh 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4GenericMessenger;

class PhysicsList : public G4VModularPhysicsList
{
 public:
  // emOption: 3 -> G4EmStandardPhysics_option3, 4 -> option4 (default)
  explicit PhysicsList(G4int emOption = 4);
  ~PhysicsList() override;

  void SetCuts() override;

  // Invoked by the /physics/absorberCut messenger
  void SetAbsorberCut(G4double value);

 private:
  G4double fAbsorberCut;  // production cut (length) for AbsorberRegion
  G4GenericMessenger* fMessenger = nullptr;
};

#endif
