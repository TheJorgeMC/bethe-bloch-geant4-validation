function dZA = betheBlochDeltaZA(material)
%BETHEBLOCHDELTAZA Absolute uncertainty on Z/A for `material`, in mol/g,
%propagated from the atomic-weight uncertainty Delta(A) (IUPAC/CIAAW).
%
%   For pure elements: Z/A = Z / A_molar, so d(Z/A)/dA_molar =
%   -Z/A_molar^2 = -(Z/A)/A_molar, giving
%   Delta(Z/A) = (Z/A)/A_molar * Delta(A_molar).
%
%   For water (a compound, material.A_g_mol is NaN): Z/A is computed
%   via the Bragg additivity rule from the H and O atomic weights, so
%   the uncertainty is propagated through BOTH constituent elements:
%
%       Z/A = w_H*(Z_H/A_H) + w_O*(Z_O/A_O)
%       Delta(Z/A) = sqrt( (w_H*Z_H/A_H^2 * Delta(A_H))^2
%                         + (w_O*Z_O/A_O^2 * Delta(A_O))^2 )
%
%   with w_H, w_O taken as exact (NIST STAR mass fractions).

    if isnan(material.A_g_mol)
        w = betheBlochWaterComposition();
        dZA_dAH = -w.w_H * w.Z_H / w.A_H^2;
        dZA_dAO = -w.w_O * w.Z_O / w.A_O^2;
        dZA = sqrt((dZA_dAH * w.Delta_A_H)^2 + (dZA_dAO * w.Delta_A_O)^2);
        return;
    end

    dZA = (material.Z_over_A / material.A_g_mol) * material.Delta_A_g_mol;

end
