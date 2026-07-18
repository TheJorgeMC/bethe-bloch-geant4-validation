function Delta_S_linear = betheBlochPropagateLinearUncertainty(T_MeV, materialKey, varargin)
%BETHEBLOCHPROPAGATELINEARUNCERTAINTY Uncertainty Delta(S_linear) on
%the LINEAR stopping power [MeV/cm].
%
%   Since S_linear = S_mass * rho, and rho is uncorrelated with the
%   mass-stopping-power sources (K, Z/A, m_e, I), the product rule gives:
%
%       Delta(S_linear) = sqrt( (rho * Delta(S_mass))^2
%                                + (S_mass * Delta(rho))^2 )
%
%   Delta(rho) is read directly from the Material struct (PubChem
%   central value, 0.5% ASSUMED uncertainty -- see betheBlochMaterials.m).

    p = inputParser;
    addParameter(p, 'z', 1);
    addParameter(p, 'particleMassMeV', []);
    parse(p, varargin{:});

    materials = betheBlochMaterials();
    material = materials.(materialKey);

    extraArgs = {'z', p.Results.z, 'particleMassMeV', p.Results.particleMassMeV};

    S_mass = betheBlochMassStoppingPower(T_MeV, material, extraArgs{:});
    Delta_S_mass = betheBlochPropagateMassUncertainty(T_MeV, materialKey, extraArgs{:});

    rho = material.rho_g_cm3;
    d_rho = material.Delta_rho_g_cm3;

    Delta_S_linear = sqrt((rho .* Delta_S_mass).^2 + (S_mass .* d_rho).^2);

end
