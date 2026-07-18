function [Delta_S, breakdown] = betheBlochPropagateMassUncertainty(T_MeV, materialKey, varargin)
%BETHEBLOCHPROPAGATEMASSUNCERTAINTY Uncertainty Delta(S) on the mass
%stopping power S(T) [MeV cm^2/g], via the standard uncorrelated
%error-propagation formula:
%
%       Delta(f) = sqrt( sum_i (df/dx_i * Delta(x_i))^2 )
%
%   with x_i in {K, Z/A, m_e, I}. Delta(I) and Delta(Z/A) (via
%   Delta(A)) are read directly from the `Material` struct -- see
%   betheBlochMaterials.m for sourcing (ICRU 90 / ICRU 49 for I,
%   IUPAC/CIAAW for A).
%
%   [Delta_S, breakdown] = betheBlochPropagateMassUncertainty(...)
%   also returns a struct `breakdown` with fields K, ZA, me, I holding
%   the individual per-source terms (same units; these combine in
%   QUADRATURE, not by direct summation, to give Delta_S).
%
%   Optional name/value pairs: 'z' (default 1), 'particleMassMeV'
%   (default proton mass).

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

    materials = betheBlochMaterials();
    material = materials.(materialKey);

    [beta, gamma] = betheBlochBetaGamma(T_MeV, particleMassMeV);

    K_val = c.K_MEV_CM2_PER_MOL;
    ZA_val = material.Z_over_A;
    me_val = c.ELECTRON_MASS_MEV;
    I_val = material.I_eV * c.EV_TO_MEV;

    dK = c.DELTA_K_MEV_CM2_PER_MOL;
    dZA = betheBlochDeltaZA(material);
    dme = c.DELTA_ELECTRON_MASS_MEV;
    dI = material.Delta_I_eV * c.EV_TO_MEV;

    d = betheBlochPartialDerivatives(K_val, ZA_val, me_val, I_val, beta, gamma, z);

    term_K = d.dK .* dK;
    term_ZA = d.dZA .* dZA;
    term_me = d.dme .* dme;
    term_I = d.dI .* dI;

    Delta_S = sqrt(term_K.^2 + term_ZA.^2 + term_me.^2 + term_I.^2);

    if nargout > 1
        breakdown.K = term_K;
        breakdown.ZA = term_ZA;
        breakdown.me = term_me;
        breakdown.I = term_I;
    end

end
