// ============================================================================
// Run.cc
// ============================================================================
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4EmCalculator.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include <cfloat>
#include <cmath>

Run::Run(const DetectorConstruction* det) : G4Run(), fDetector(det) {}

void Run::SetPrimary(const G4ParticleDefinition* particle, G4double energy)
{
  fParticle = particle;
  fEkin = energy;
}

void Run::AddEvent(G4double edepPrimary, G4double edepSecondary,
                   G4double escapedSecondary, G4double exitEnergy,
                   G4double trackLength)
{
  const G4double edepTotal = edepPrimary + edepSecondary;
  fSumEdepTotal += edepTotal;
  fSumEdepTotal2 += edepTotal * edepTotal;
  fSumEdepPrimary += edepPrimary;
  fSumEdepSecondary += edepSecondary;
  fSumEscapedSecondary += escapedSecondary;
  if (exitEnergy >= 0.) {
    fSumExitEnergy += exitEnergy;
    ++fNExit;
  }
  fSumTrackLength += trackLength;
  ++fNEvents;
}

void Run::Merge(const G4Run* aRun)
{
  // Fusion thread-safe: G4 la invoca secuencialmente sobre el Run del master
  // con cada Run local de los hilos de trabajo.
  const Run* localRun = static_cast<const Run*>(aRun);

  if (fParticle == nullptr) {
    fParticle = localRun->fParticle;
    fEkin = localRun->fEkin;
  }
  fSumEdepTotal += localRun->fSumEdepTotal;
  fSumEdepTotal2 += localRun->fSumEdepTotal2;
  fSumEdepPrimary += localRun->fSumEdepPrimary;
  fSumEdepSecondary += localRun->fSumEdepSecondary;
  fSumEscapedSecondary += localRun->fSumEscapedSecondary;
  fSumExitEnergy += localRun->fSumExitEnergy;
  fNExit += localRun->fNExit;
  fSumTrackLength += localRun->fSumTrackLength;
  fNEvents += localRun->fNEvents;

  G4Run::Merge(aRun);
}

void Run::EndOfRun() const
{
  if (fNEvents == 0 || fParticle == nullptr) {
    G4cout << "### Run::EndOfRun: sin eventos, nada que reportar" << G4endl;
    return;
  }

  G4Material* material = fDetector->GetAbsorberMaterial();
  const G4double thickness = fDetector->GetThickness();
  const G4double n = static_cast<G4double>(fNEvents);

  // --- Estadisticas de la simulacion -----------------------------------
  const G4double meanEdep = fSumEdepTotal / n;
  G4double rms2 = fSumEdepTotal2 / n - meanEdep * meanEdep;
  const G4double rmsEdep = (rms2 > 0.) ? std::sqrt(rms2) : 0.;
  const G4double meanEdepPrim = fSumEdepPrimary / n;
  const G4double meanEdepSec = fSumEdepSecondary / n;
  const G4double meanEscaped = fSumEscapedSecondary / n;
  const G4double meanTrackLen = fSumTrackLength / n;
  const G4double meanExitAll = fSumExitEnergy / n;  // 0 para los que no salen
  const G4double meanExit =
      (fNExit > 0) ? fSumExitEnergy / static_cast<G4double>(fNExit) : 0.;

  // dE/dx "local" simulado: energia depositada por el PRIMARIO por unidad de
  // espesor. Aproxima el dE/dx RESTRINGIDO (la perdida continua del proceso
  // de ionizacion) mas los depositos discretos que quedan bajo umbral.
  const G4double dedxPrimSim = meanEdepPrim / thickness;

  // dE/dx "total" simulado: perdida de energia del primario por unidad de
  // espesor, incluyendo la energia cedida a los rayos delta explicitos
  // (= Einc - Esalida, promediada). Aproxima el dE/dx NO RESTRINGIDO cuando
  // el slab es delgado.
  const G4double dedxTotalSim =
      (fNExit == fNEvents) ? (fEkin - meanExit) / thickness : 0.;

  // --- Valores de referencia de G4EmCalculator (patron TestEm0) ---------
  G4EmCalculator emCal;
  const G4Region* region = G4RegionStore::GetInstance()->GetRegion(
      DetectorConstruction::kRegionName, false);

  // Restringido: dE/dx tabulado del proceso de ionizacion con el corte
  // activo en la region del slab.
  const G4double dedxRestricted =
      emCal.GetDEDX(fEkin, fParticle, material, region);

  // No restringido: dE/dx electronico calculado con T_cut = DBL_MAX
  // (es decir, integrando hasta T_max). Esta es la cantidad comparable con
  // la formula de Bethe-Bloch "de libro de texto".
  const G4double dedxUnrestricted =
      emCal.ComputeElectronicDEDX(fEkin, fParticle, material, DBL_MAX);

  // Rango CSDA (requiere G4EmParameters::SetBuildCSDARange(true), activado
  // en el constructor de PhysicsList).
  const G4double csdaRange = emCal.GetCSDARange(fEkin, fParticle, material);

  // --- Balance de energia por evento ------------------------------------
  // Einc ≈ Edep(prim) + Edep(sec) + E(escapada por sec) + E(salida prim).
  // El residuo refleja canales no contabilizados (p.ej. fotones de
  // fluorescencia que escapan, redondeo numerico) y debe ser pequeño.
  const G4double balance =
      meanEdepPrim + meanEdepSec + meanEscaped + meanExitAll;
  const G4double residual = fEkin - balance;

  G4cout << "\n============================================================\n"
         << " Resumen del run — validacion Bethe-Bloch\n"
         << "============================================================\n"
         << " Primario           : " << fParticle->GetParticleName() << ", E = "
         << G4BestUnit(fEkin, "Energy") << "\n"
         << " Material           : " << material->GetName() << "\n"
         << " Espesor del slab   : " << G4BestUnit(thickness, "Length") << "\n"
         << " Eventos            : " << fNEvents << "\n"
         << "------------------------------------------------------------\n"
         << " <Edep total>       : " << G4BestUnit(meanEdep, "Energy")
         << "  (rms/straggling: " << G4BestUnit(rmsEdep, "Energy") << ")\n"
         << " <Edep primario>    : " << G4BestUnit(meanEdepPrim, "Energy") << "\n"
         << " <Edep secundarios> : " << G4BestUnit(meanEdepSec, "Energy") << "\n"
         << " <E escapada (sec)> : " << G4BestUnit(meanEscaped, "Energy") << "\n"
         << " Primarios salientes: " << fNExit << " / " << fNEvents;
  if (fNExit > 0) {
    G4cout << "   <E salida> = " << G4BestUnit(meanExit, "Energy");
  }
  G4cout << "\n <longitud de traza>: " << G4BestUnit(meanTrackLen, "Length")
         << "\n"
         << "------------------------------------------------------------\n"
         << " dE/dx sim (local, ~restringido)   : "
         << G4BestUnit(dedxPrimSim, "Energy/Length") << "\n";
  if (dedxTotalSim > 0.) {
    G4cout << " dE/dx sim (Einc-Esal, ~total)     : "
           << G4BestUnit(dedxTotalSim, "Energy/Length") << "\n";
  } else {
    G4cout << " dE/dx sim (Einc-Esal, ~total)     : n/a (no todos los "
              "primarios salieron del slab)\n";
  }
  G4cout << " G4EmCalculator dE/dx restringido  : "
         << G4BestUnit(dedxRestricted, "Energy/Length") << "\n"
         << " G4EmCalculator dE/dx NO restringido: "
         << G4BestUnit(dedxUnrestricted, "Energy/Length") << "\n"
         << " G4EmCalculator rango CSDA         : "
         << G4BestUnit(csdaRange, "Length") << "\n"
         << "------------------------------------------------------------\n"
         << " Balance <por evento>: Edep_p + Edep_s + E_esc + E_sal = "
         << G4BestUnit(balance, "Energy") << "\n"
         << " Residuo (Einc - balance)          : "
         << G4BestUnit(residual, "Energy") << "  ("
         << 100. * residual / fEkin << " % de Einc)\n"
         << "============================================================\n"
         << G4endl;
}
