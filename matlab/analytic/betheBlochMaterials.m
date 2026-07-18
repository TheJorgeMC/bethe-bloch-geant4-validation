function materials = betheBlochMaterials()
%BETHEBLOCHMATERIALS Properties of the 4 absorber materials used
%throughout the project, with I, Delta_I, rho, Delta_rho, A and
%Delta_A as EXPLICIT fields (mirrors python/analytic Material
%dataclass exactly -- no source is treated as "negligible" and
%hidden away; everything is a first-class, cited field).
%
%   ------------------------------------------------------------------
%   Mean excitation energies I and their uncertainties:
%     Water:     ICRU Report 90 (2016), Sec. 5.3 "Mean Excitation Energy
%                of Liquid Water" -- recommended value 78(2) eV.
%     Al/Cu/Pb:  ICRU Report 49 (1993), Table 2.8, Sec. 2.5.1 "Mean
%                Excitation Energies for Elements".
%
%   Standard atomic weights A and their uncertainties:
%     IUPAC/CIAAW, Table 2 "List of Elements in Atomic Number Order"
%     (https://iupac.qmul.ac.uk/AtWt/), 2021 recommended values.
%     NOTE: Pb's standard atomic weight is 207.2(1.1) g/mol -- the
%     parenthetical uncertainty is +/-1.1, NOT +/-0.1. This reflects
%     genuine natural isotopic variability (Pb's standard atomic weight
%     was changed to an interval, [206.14, 207.94], in the 2021 IUPAC
%     revision); 207.2(1.1) is IUPAC's own conservative symmetric
%     value chosen to cover that full interval.
%
%   Densities and their uncertainties:
%     PubChem Periodic Table (https://pubchem.ncbi.nlm.nih.gov/ptable/density/).
%     PubChem does not publish a density uncertainty; a representative
%     0.5% "engineering-grade" tolerance is ASSUMED for all four
%     materials (no uncertainty source was located for solid/liquid
%     densities in the literature consulted).
%   ------------------------------------------------------------------
%
%   Returns a struct with one field per material; each field has:
%       .name              Display name
%       .Z_over_A          Effective Z/A [mol/g]
%       .I_eV              Mean excitation energy [eV]
%       .Delta_I_eV        Uncertainty on I [eV]
%       .rho_g_cm3         Density [g/cm^3]
%       .Delta_rho_g_cm3   Uncertainty on density [g/cm^3]
%       .A_g_mol           Standard atomic weight [g/mol] (NaN for water)
%       .Delta_A_g_mol     Uncertainty on A [g/mol] (NaN for water)
%       .matno             NIST STAR material number (traceability)

    materials.aluminium = makeMaterial('ALUMINIUM', 13.0/26.9815384, ...
        166.0, 2.0, ...                     % ICRU 49, Table 2.8
        2.70, 2.70*0.005, ...               % PubChem; 0.5% ASSUMED
        26.9815384, 0.0000003, ...          % IUPAC: 26.981 5384(3)
        '013');

    materials.copper = makeMaterial('COPPER', 29.0/63.546, ...
        322.0, 10.0, ...                    % ICRU 49, Table 2.8
        8.96, 8.96*0.005, ...               % PubChem; 0.5% ASSUMED
        63.546, 0.003, ...                  % IUPAC: 63.546(3)
        '029');

    % NOTE on the Pb atomic weight uncertainty below (1.1, not 0.1):
    % see the header comment above for the full explanation.
    materials.lead = makeMaterial('LEAD', 82.0/207.2, ...
        823.0, 30.0, ...                    % ICRU 49, Table 2.8
        11.34, 11.34*0.005, ...             % PubChem; 0.5% ASSUMED
        207.2, 1.1, ...                     % IUPAC: 207.2(1.1)
        '082');

    % Liquid water: compound, Z/A by the Bragg rule (see
    % betheBlochWaterComposition.m for the constituent H/O values).
    w = betheBlochWaterComposition();
    water_Z_over_A = w.w_H * (w.Z_H / w.A_H) + w.w_O * (w.Z_O / w.A_O);

    materials.water = makeMaterial('WATER', water_Z_over_A, ...
        78.0, 2.0, ...                      % ICRU 90 (2016), Sec. 5.3
        1.0, 1.0*0.005, ...                 % PubChem; 0.5% ASSUMED
        NaN, NaN, ...                       % compound -- see betheBlochWaterComposition.m
        '276');

end

function m = makeMaterial(name, Z_over_A, I_eV, Delta_I_eV, ...
        rho_g_cm3, Delta_rho_g_cm3, A_g_mol, Delta_A_g_mol, matno)
    m.name = name;
    m.Z_over_A = Z_over_A;
    m.I_eV = I_eV;
    m.Delta_I_eV = Delta_I_eV;
    m.rho_g_cm3 = rho_g_cm3;
    m.Delta_rho_g_cm3 = Delta_rho_g_cm3;
    m.A_g_mol = A_g_mol;
    m.Delta_A_g_mol = Delta_A_g_mol;
    m.matno = matno;
end