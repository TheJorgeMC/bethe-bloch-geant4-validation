function [beta, gamma] = betheBlochBetaGamma(T_MeV, particleMassMeV)
%BETHEBLOCHBETAGAMMA Relativistic beta and gamma from kinetic energy
%T_MeV (MeV) of a particle with rest mass energy particleMassMeV (MeV).
%
%   gamma = 1 + T / (M c^2)
%   beta  = sqrt(1 - 1/gamma^2)
%
%   If particleMassMeV is empty or omitted, the proton mass is used.

    if nargin < 2 || isempty(particleMassMeV)
        c = betheBlochConstants();
        particleMassMeV = c.PROTON_MASS_MEV;
    end

    T_MeV = T_MeV(:)';
    gamma = 1.0 + T_MeV ./ particleMassMeV;
    beta = sqrt(1.0 - 1.0 ./ gamma.^2);

end
