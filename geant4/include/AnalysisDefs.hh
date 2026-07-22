// ============================================================================
// AnalysisDefs.hh
// Identificadores compartidos de histogramas y columnas del ntuple, para que
// RunAction (que los crea) y EventAction (que los llena) no dependan de
// literales duplicados.
// ============================================================================
#ifndef AnalysisDefs_hh
#define AnalysisDefs_hh 1

#include "globals.hh"

namespace analysis
{
// --- Histogramas 1D ---
// H1 0: energia depositada vs profundidad z en el slab (curva de Bragg si
//       /absorber/numberOfLayers > 1). Se rebinnea en cada BeginOfRunAction.
constexpr G4int kH1Depth = 0;

// --- Columnas del ntuple "slab" (una fila por evento, todas en MeV/mm) ---
enum NtupleColumns : G4int
{
  kColEIncident = 0,       // energia cinetica del proton primario al nacer
  kColEdepPrimary = 1,     // energia depositada en el slab por el primario
  kColEdepSecondary = 2,   // energia depositada en el slab por secundarios
  kColEEscapedSec = 3,     // energia cinetica sacada del slab por secundarios
  kColEExitPrimary = 4,    // energia cinetica del primario al salir (-1 si no sale)
  kColTrackLenPrimary = 5  // longitud de traza del primario dentro del slab
};
}  // namespace analysis

#endif
