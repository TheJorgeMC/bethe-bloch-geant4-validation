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
  // Thread-safe merge: G4 invokes this sequentially on the master's Run with
  // each worker thread's local Run.
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
    G4cout << "### Run::EndOfRun: no events, nothing to report" << G4endl;
    return;
  }

  G4Material* material = fDetector->GetAbsorberMaterial();
  const G4double thickness = fDetector->GetThickness();
  const G4double n = static_cast<G4double>(fNEvents);

  // --- Simulation statistics ---------------------------------------------
  const G4double meanEdep = fSumEdepTotal / n;
  G4double rms2 = fSumEdepTotal2 / n - meanEdep * meanEdep;
  const G4double rmsEdep = (rms2 > 0.) ? std::sqrt(rms2) : 0.;
  const G4double meanEdepPrim = fSumEdepPrimary / n;
  const G4double meanEdepSec = fSumEdepSecondary / n;
  const G4double meanEscaped = fSumEscapedSecondary / n;
  const G4double meanTrackLen = fSumTrackLength / n;
  const G4double meanExitAll = fSumExitEnergy / n;  // 0 for non-exiting ones
  const G4double meanExit =
      (fNExit > 0) ? fSumExitEnergy / static_cast<G4double>(fNExit) : 0.;

  // Simulated "local" dE/dx: energy deposited by the PRIMARY per unit
  // thickness. Approximates the RESTRICTED dE/dx (the continuous energy loss
  // of the ionisation process) plus below-threshold discrete deposits.
  const G4double dedxPrimSim = meanEdepPrim / thickness;

  // Simulated "total" dE/dx: primary energy loss per unit thickness,
  // including the energy handed to explicit delta rays
  // (= Einc - Eexit, averaged). Approximates the UNRESTRICTED dE/dx for a
  // thin slab.
  const G4double dedxTotalSim =
      (fNExit == fNEvents) ? (fEkin - meanExit) / thickness : 0.;

  // --- Reference values from G4EmCalculator (TestEm0 pattern) ------------
  G4EmCalculator emCal;
  const G4Region* region = G4RegionStore::GetInstance()->GetRegion(
      DetectorConstruction::kRegionName, false);

  // Restricted: tabulated dE/dx of the ionisation process with the cut
  // active in the slab region.
  const G4double dedxRestricted =
      emCal.GetDEDX(fEkin, fParticle, material, region);

  // Unrestricted: electronic dE/dx computed with T_cut = DBL_MAX (i.e.
  // integrating up to T_max). This is the quantity comparable with the
  // "textbook" Bethe-Bloch formula.
  const G4double dedxUnrestricted =
      emCal.ComputeElectronicDEDX(fEkin, fParticle, material, DBL_MAX);

  // CSDA range (requires G4EmParameters::SetBuildCSDARange(true), enabled
  // in the PhysicsList constructor).
  const G4double csdaRange = emCal.GetCSDARange(fEkin, fParticle, material);

  // --- Per-event energy balance -------------------------------------------
  // Einc ~ Edep(prim) + Edep(sec) + E(escaped via sec) + E(primary exit).
  // The residual reflects unaccounted channels (e.g. fluorescence photons
  // escaping, numerical rounding) and must be small.
  const G4double balance =
      meanEdepPrim + meanEdepSec + meanEscaped + meanExitAll;
  const G4double residual = fEkin - balance;

  G4cout << "\n============================================================\n"
         << " Run summary — Bethe-Bloch validation\n"
         << "============================================================\n"
         << " Primary            : " << fParticle->GetParticleName() << ", E = "
         << G4BestUnit(fEkin, "Energy") << "\n"
         << " Material           : " << material->GetName() << "\n"
         << " Slab thickness     : " << G4BestUnit(thickness, "Length") << "\n"
         << " Events             : " << fNEvents << "\n"
         << "------------------------------------------------------------\n"
         << " <Edep total>       : " << G4BestUnit(meanEdep, "Energy")
         << "  (rms/straggling: " << G4BestUnit(rmsEdep, "Energy") << ")\n"
         << " <Edep primary>     : " << G4BestUnit(meanEdepPrim, "Energy") << "\n"
         << " <Edep secondaries> : " << G4BestUnit(meanEdepSec, "Energy") << "\n"
         << " <E escaped (sec)>  : " << G4BestUnit(meanEscaped, "Energy") << "\n"
         << " Exiting primaries  : " << fNExit << " / " << fNEvents;
  if (fNExit > 0) {
    G4cout << "   <E exit> = " << G4BestUnit(meanExit, "Energy");
  }
  G4cout << "\n <track length>     : " << G4BestUnit(meanTrackLen, "Length")
         << "\n"
         << "------------------------------------------------------------\n"
         << " dE/dx sim (local, ~restricted)     : "
         << G4BestUnit(dedxPrimSim, "Energy/Length") << "\n";
  if (dedxTotalSim > 0.) {
    G4cout << " dE/dx sim (Einc-Eexit, ~total)     : "
           << G4BestUnit(dedxTotalSim, "Energy/Length") << "\n";
  } else {
    G4cout << " dE/dx sim (Einc-Eexit, ~total)     : n/a (not all primaries "
              "exited the slab)\n";
  }
  G4cout << " G4EmCalculator restricted dE/dx    : "
         << G4BestUnit(dedxRestricted, "Energy/Length") << "\n"
         << " G4EmCalculator UNRESTRICTED dE/dx  : "
         << G4BestUnit(dedxUnrestricted, "Energy/Length") << "\n"
         << " G4EmCalculator CSDA range          : "
         << G4BestUnit(csdaRange, "Length") << "\n"
         << "------------------------------------------------------------\n"
         << " Balance <per event>: Edep_p + Edep_s + E_esc + E_exit = "
         << G4BestUnit(balance, "Energy") << "\n"
         << " Residual (Einc - balance)          : "
         << G4BestUnit(residual, "Energy") << "  ("
         << 100. * residual / fEkin << " % of Einc)\n"
         << "============================================================\n"
         << G4endl;
}
