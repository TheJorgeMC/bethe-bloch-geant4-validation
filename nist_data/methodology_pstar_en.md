# NIST PSTAR Methodology — Reference Note

**Source:** NIST Standard Reference Database 124, PSTAR
(`https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html`)

PSTAR computes proton stopping power and range by combining two energy
regimes, each treated with a different approach:

1. **Low energy (roughly T ≲ 1 MeV, precisely the region where Bethe
   theory loses validity):** the electronic stopping power is obtained by
   **cubic spline interpolation** over experimental data compiled and
   evaluated by ICRU (Report 49) and related compilations (Andersen &
   Ziegler, Ziegler). In this regime, Born-Bethe theory is unreliable
   because the proton velocity is comparable to the orbital velocities of
   the target's atomic electrons, so semi-empirical fits to measured data
   are used instead.

2. **High energy (above the range reliably covered by the empirical
   fits):** the **Bethe formula** is used (the same functional form
   derived from first-order Born theory, with the relativistic logarithm
   and the 1/β² prefactor), including the standard corrections from the
   ICRU 37 / ICRU 90 formalism:
   - **Shell correction** — relevant in the regime where the Born
     approximation starts to break down due to the absorber's atomic
     structure (important for Pb across nearly the whole clinical range,
     as already established in the earlier derivation work).
   - **Barkas-Andersen correction** (z³ term, charge asymmetry).
   - **Bloch correction** (interpolation between the Born limit and the
     classical Bohr limit, governed by the Sommerfeld parameter
     η = zα/β).

   Nuclear stopping power (STOP(n)) is calculated separately using an
   elastic nucleus-nucleus collision model (universal screened
   potential), dominant only at very low energies (keV range) and
   negligible in the 3–300 MeV clinical range.

3. **Matching:** the two regions are joined continuously in both energy
   and derivative at the crossover point, so the final table produced by
   PSTAR is a smooth function with no discontinuities, even though it
   internally combines a spline fit (low T) with a corrected analytic
   formula (high T).

**Relevance to this project:** across the clinical proton therapy range
(3–300 MeV, βγ ≳ 0.08), the proton is firmly within the regime dominated
by the Bethe formula with ICRU 90 corrections, making PSTAR a suitable
and theoretically consistent reference for the Bethe-Bloch formulation
(relativistic Bethe-Bloch + shell + Barkas/Bloch corrections) being
validated in this work — as opposed to the very-low-energy experimental
spline regime, which falls outside the range of interest here.

**Suggested citation:**
Berger, M.J., Coursey, J.S., Zucker, M.A., and Chang, J. (2005),
*Stopping-Power and Range Tables for Electrons, Protons, and Helium Ions*
(NIST Standard Reference Database 124), National Institute of Standards
and Technology, Gaithersburg, MD. https://dx.doi.org/10.18434/T4NC7P
