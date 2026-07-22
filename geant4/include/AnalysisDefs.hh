// ============================================================================
// AnalysisDefs.hh
// Shared histogram and ntuple-column identifiers, so that RunAction (which
// creates them) and EventAction (which fills them) do not rely on duplicated
// literals.
// ============================================================================
#ifndef AnalysisDefs_hh
#define AnalysisDefs_hh 1

#include "globals.hh"

namespace analysis
{
// --- 1D histograms ---
// H1 0: energy deposited vs depth z in the slab (Bragg curve when
//       /absorber/numberOfLayers > 1). Rebinned at every BeginOfRunAction.
constexpr G4int kH1Depth = 0;

// --- Columns of the "slab" ntuple (one row per event, all in MeV/mm) ---
enum NtupleColumns : G4int
{
  kColEIncident = 0,       // primary proton kinetic energy at birth
  kColEdepPrimary = 1,     // energy deposited in the slab by the primary
  kColEdepSecondary = 2,   // energy deposited in the slab by secondaries
  kColEEscapedSec = 3,     // kinetic energy carried out of the slab by secondaries
  kColEExitPrimary = 4,    // primary kinetic energy on exit (-1 if it never exits)
  kColTrackLenPrimary = 5  // primary track length inside the slab
};
}  // namespace analysis

#endif
