function c = betheBlochConstants()
%BETHEBLOCHCONSTANTS Physical constants for the analytic Bethe-Bloch
%implementation (relativistic, no corrections).
%
%   Mirrors python/analytic constants exactly, so cross-validation
%   matches to floating-point precision. Source: CODATA 2022 / PDG 2024.

    %   Electron rest mass energy m_e*c^2 in MeV (CODATA, 2022)
    c.ELECTRON_MASS_MEV = 0.51099895069;

    %   Proton rest mass energy M_p*c^2 in MeV (CODATA, 2022)
    c.PROTON_MASS_MEV = 938.27208943;

    %   PDG's K constant: K = 4*pi*N_A*r_e^2*m_e*c^2, in MeV*cm^2/mol
    %   (previously verified against own theoretical derivation, Sect. 2.6 of
    %   Bethe_Full_Derivation.pdf file in the report directory, coinciding
    %   with PDG 2024 Sect 34.2.1)
    c.K_MEV_CM2_PER_MOL = 0.307075;

    %   Conversion factor eV -> MeV
    c.EV_TO_MEV = 1e-6;

    %   Delta(m_e c^2): CODATA 2022 (NIST, https://physics.nist.gov/cgi-bin/cuu/Value?mec2mev),
    %   m_e c^2 = 0.510 998 950 69(16) MeV -- i.e. Delta = 1.6e-10 MeV
    %   (relative uncertainty ~3.1e-10).
    c.DELTA_ELECTRON_MASS_MEV = 1.6e-10;

    %   Delta(K): K = 4*pi*N_A*r_e^2*m_e*c^2 is not independently tabulated
    %   by CODATA or PDG; it is PROPAGATED from the same set of CODATA 2022
    %   constants that define it. N_A is EXACT (2019 SI redefinition); m_e c^2
    %   is known to ~3.1e-10 relative (see above) and r_e (classical electron
    %   radius) to ~4.6e-10 relative, giving K a relative uncertainty of order
    %   1e-9 -- negligible next to Delta(I), but propagated explicitly rather
    %   than assumed zero.
    c.DELTA_K_MEV_CM2_PER_MOL = 3.0e-10;

end
