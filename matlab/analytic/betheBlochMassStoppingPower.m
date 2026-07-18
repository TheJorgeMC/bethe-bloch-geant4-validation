function S = betheBlochMassStoppingPower(T_MeV, material, varargin)
%BETHEBLOCHMASSSTOPPINGPOWER Bethe-Bloch relativistic mass stopping
%power, SIN shell/Barkas/Bloch/density corrections.
%
%   S = betheBlochMassStoppingPower(T_MeV, material) returns
%   -(1/rho)(dE/dx), in MeV*cm^2/g, for a proton (z=1) with kinetic
%   energy T_MeV (MeV, scalar or vector) in `material` (struct field
%   from betheBlochMaterials.m).
%
%   S = betheBlochMassStoppingPower(T_MeV, material, 'z', z, ...
%       'particleMassMeV', M) sets projectile charge/mass (default:
%       z=1, M = proton mass).
%
%       S(T) = K * z^2 * (Z/A) * (1/beta^2) *
%              [ ln(2*m_e*c^2*beta^2*gamma^2 / I) - beta^2 ]
%
%   ALCANCE: does NOT include shell, Barkas, Bloch, or density-effect
%   corrections. Expected to underestimate NIST PSTAR / Geant4,
%   increasingly so at low energy and high Z. See the theoretical
%   derivation document, Sec. "Project Scope and Future Work".

    p = inputParser;
    addParameter(p, 'z', 1);
    addParameter(p, 'particleMassMeV', []);
    parse(p, varargin{:});
    z = p.Results.z;

    c = betheBlochConstants();
    if isempty(p.Results.particleMassMeV)
        particleMassMeV = c.PROTON_MASS_MEV;
    else
        particleMassMeV = p.Results.particleMassMeV;
    end

    [beta, gamma] = betheBlochBetaGamma(T_MeV, particleMassMeV);

    I_MeV = material.I_eV * c.EV_TO_MEV;

    log_arg = (2.0 .* c.ELECTRON_MASS_MEV .* beta.^2 .* gamma.^2) ./ I_MeV;
    main_bracket = log(log_arg) - beta.^2;

    S = c.K_MEV_CM2_PER_MOL .* z^2 .* material.Z_over_A .* main_bracket ./ beta.^2;

end
