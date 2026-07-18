function S_linear = betheBlochLinearStoppingPower(T_MeV, material, varargin)
%BETHEBLOCHLINEARSTOPPINGPOWER Linear stopping power -(dE/dx), in
%MeV/cm, for a particle with kinetic energy T_MeV in `material`. This
%is simply the mass stopping power times the material density.
%
%   See betheBlochMassStoppingPower.m for optional 'z'/'particleMassMeV'.

    S_mass = betheBlochMassStoppingPower(T_MeV, material, varargin{:});
    S_linear = S_mass .* material.rho_g_cm3;

end
